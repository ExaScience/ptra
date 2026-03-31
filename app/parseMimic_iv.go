// PTRA: Patient Trajectory Analysis Library
// Copyright (c) 2022 imec vzw.

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as
// published by the Free Software Foundation, either version 3 of the
// License, or (at your option) any later version, and Additional Terms
// (see below).

// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Affero General Public License for more details.

// You should have received a copy of the GNU Affero General Public
// License and Additional Terms along with this program. If not, see
// <https://github.com/ExaScience/ptra/blob/master/LICENSE.txt>.

package app

import (
	"encoding/csv"
	"encoding/json"
	"fmt"
	"io"
	"io/ioutil"
	"math"
	"os"
	"path/filepath"
	"ptra/trajectory"
	"ptra/utils"
	"sort"
	"strconv"
	"strings"
)

// ---------------------------------------------------------------------------
// Shared MIMIC-IV helpers
// ---------------------------------------------------------------------------

// ParseMIMIC4Date parses a date string in YYYY-MM-DD format into a DiagnosisDate.
func ParseMIMIC4Date(dateStr string) (trajectory.DiagnosisDate, bool) {
	parts := strings.Split(dateStr, "-")
	if len(parts) != 3 {
		return trajectory.DiagnosisDate{}, false
	}
	year, err := strconv.Atoi(parts[0])
	if err != nil {
		return trajectory.DiagnosisDate{}, false
	}
	month, err := strconv.Atoi(parts[1])
	if err != nil {
		return trajectory.DiagnosisDate{}, false
	}
	day, err := strconv.Atoi(parts[2])
	if err != nil {
		return trajectory.DiagnosisDate{}, false
	}
	return trajectory.DiagnosisDate{Year: year, Month: month, Day: day}, true
}

// PatientDemo holds demographic data for a MIMIC-IV patient.
type PatientDemo struct {
	Sex int
	Yob int
}

// LoadMIMIC4PatientDemographics loads patients.csv and returns a map subject_id -> PatientDemo.
func LoadMIMIC4PatientDemographics(patientsFile string) map[string]*PatientDemo {
	fmt.Println("Parsing MIMIC-IV patient demographics from: ", patientsFile)
	pFile, err := os.Open(patientsFile)
	if err != nil {
		panic(err)
	}
	pReader := csv.NewReader(pFile)
	if _, err := pReader.Read(); err != nil {
		panic(err)
	}
	demoMap := map[string]*PatientDemo{}
	for {
		record, err := pReader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			panic(err)
		}
		subjectIDStr := record[0]
		sex := trajectory.Male
		if record[1] == "F" {
			sex = trajectory.Female
		}
		anchorAge, _ := strconv.Atoi(record[2])
		anchorYear, _ := strconv.Atoi(record[3])
		yob := 0
		if anchorYear > 0 && anchorAge > 0 {
			yob = anchorYear - anchorAge
		}
		demoMap[subjectIDStr] = &PatientDemo{Sex: sex, Yob: yob}
	}
	pFile.Close()
	fmt.Println("Loaded demographics for ", len(demoMap), " patients.")
	return demoMap
}

// ---------------------------------------------------------------------------
// HCPCS parser
// ---------------------------------------------------------------------------

// ParseHCPCSEventsWithPatients parses patients.csv then streams hcpcsevents.csv
// row-by-row to build a trajectory.Experiment.
func ParseHCPCSEventsWithPatients(hcpcsFile, patientsFile string, nofCohortAges int,
	filters []trajectory.PatientFilter) (*trajectory.Experiment, *trajectory.PatientMap) {

	demoMap := LoadMIMIC4PatientDemographics(patientsFile)

	fmt.Println("Parsing MIMIC-IV HCPCS events from: ", hcpcsFile)
	file, err := os.Open(hcpcsFile)
	if err != nil {
		panic(err)
	}
	defer func() {
		if err := file.Close(); err != nil {
			panic(err)
		}
	}()

	reader := csv.NewReader(file)
	header, err := reader.Read()
	if err != nil {
		panic(fmt.Sprint("failed to read CSV header: ", err))
	}
	colSubjectID, colChartDate, colHcpcsCd, colShortDesc := -1, -1, -1, -1
	for i, h := range header {
		switch strings.TrimSpace(strings.ToLower(h)) {
		case "subject_id":
			colSubjectID = i
		case "chartdate":
			colChartDate = i
		case "hcpcs_cd":
			colHcpcsCd = i
		case "short_description":
			colShortDesc = i
		}
	}
	if colSubjectID == -1 || colChartDate == -1 || colHcpcsCd == -1 {
		panic("hcpcsevents.csv must contain subject_id, chartdate, and hcpcs_cd columns")
	}

	patientMap := &trajectory.PatientMap{
		PIDMap:       map[int]*trajectory.Patient{},
		PIDStringMap: map[string]int{},
	}
	hcpcsCodeToDID := map[string]int{}
	didToName := map[int]string{}
	didToCode := map[int]string{}
	didCounter := 0
	maxYOB := 1850
	minYOB := 2100
	eventCtr := 0
	skippedCtr := 0

	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			panic(err)
		}
		eventCtr++
		subjectIDStr := record[colSubjectID]
		chartDateStr := record[colChartDate]
		hcpcsCd := record[colHcpcsCd]

		date, ok := ParseMIMIC4Date(chartDateStr)
		if !ok {
			skippedCtr++
			continue
		}
		if _, exists := hcpcsCodeToDID[hcpcsCd]; !exists {
			hcpcsCodeToDID[hcpcsCd] = didCounter
			name := hcpcsCd
			if colShortDesc != -1 && colShortDesc < len(record) && record[colShortDesc] != "" {
				name = record[colShortDesc]
			}
			didToName[didCounter] = name
			didToCode[didCounter] = hcpcsCd
			didCounter++
		}
		if _, exists := patientMap.PIDStringMap[subjectIDStr]; !exists {
			patientMap.Ctr++
			pid := patientMap.Ctr
			sex := trajectory.Male
			yob := 0
			if demo, found := demoMap[subjectIDStr]; found {
				sex = demo.Sex
				yob = demo.Yob
			}
			if sex == trajectory.Male {
				patientMap.MaleCtr++
			} else {
				patientMap.FemaleCtr++
			}
			if yob > 0 {
				maxYOB = utils.MaxInt(yob, maxYOB)
				minYOB = utils.MinInt(yob, minYOB)
			}
			patient := &trajectory.Patient{
				PID: pid, PIDString: subjectIDStr, YOB: yob,
				CohortAge: 0, Sex: sex,
				Diagnoses: []*trajectory.Diagnosis{}, Region: 0,
			}
			patientMap.PIDMap[pid] = patient
			patientMap.PIDStringMap[subjectIDStr] = pid
		}
		pid := patientMap.PIDStringMap[subjectIDStr]
		patient := patientMap.PIDMap[pid]
		did := hcpcsCodeToDID[hcpcsCd]
		diagnosis := &trajectory.Diagnosis{PID: pid, DID: did, Date: date}
		trajectory.AddDiagnosis(patient, diagnosis)
	}

	fmt.Println("Parsed ", eventCtr, " HCPCS events for ", len(patientMap.PIDMap), " patients.")
	fmt.Printf("Males: %d, Females: %d\n", patientMap.MaleCtr, patientMap.FemaleCtr)
	fmt.Println("Unique HCPCS codes: ", didCounter)
	fmt.Println("Skipped events (bad date): ", skippedCtr)
	fmt.Printf("Year of birth range: %d - %d\n", minYOB, maxYOB)

	return BuildMIMIC4Experiment(patientMap, didToName, didToCode, didCounter,
		nofCohortAges, maxYOB, minYOB, "mimic4_hcpcs", filters)
}

// ---------------------------------------------------------------------------
// diagnoses_icd support
// ---------------------------------------------------------------------------

// LoadAdmissionDates loads admissions.csv and returns a map hadm_id -> DiagnosisDate.
func LoadAdmissionDates(admissionsFile string) map[string]trajectory.DiagnosisDate {
	fmt.Println("Loading admission dates from: ", admissionsFile)
	file, err := os.Open(admissionsFile)
	if err != nil {
		panic(err)
	}
	defer file.Close()
	reader := csv.NewReader(file)
	header, err := reader.Read()
	if err != nil {
		panic(err)
	}
	colHadmID, colAdmitTime := -1, -1
	for i, h := range header {
		switch strings.TrimSpace(strings.ToLower(h)) {
		case "hadm_id":
			colHadmID = i
		case "admittime":
			colAdmitTime = i
		}
	}
	if colHadmID == -1 || colAdmitTime == -1 {
		panic("admissions.csv must contain hadm_id and admittime columns")
	}
	result := map[string]trajectory.DiagnosisDate{}
	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			panic(err)
		}
		if len(record[colAdmitTime]) >= 10 {
			if date, ok := ParseMIMIC4Date(record[colAdmitTime][:10]); ok {
				result[record[colHadmID]] = date
			}
		}
	}
	fmt.Println("Loaded admission dates for ", len(result), " admissions.")
	return result
}

// LoadICDNames loads d_icd_diagnoses.csv and returns a map icd_code -> long_title.
func LoadICDNames(file string) map[string]string {
	fmt.Println("Loading ICD diagnosis names from: ", file)
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	reader := csv.NewReader(f)
	header, err := reader.Read()
	if err != nil {
		panic(err)
	}
	colCode, colTitle := -1, -1
	for i, h := range header {
		switch strings.TrimSpace(strings.ToLower(h)) {
		case "icd_code":
			colCode = i
		case "long_title":
			colTitle = i
		}
	}
	if colCode == -1 || colTitle == -1 {
		panic("d_icd_diagnoses.csv must contain icd_code and long_title columns")
	}
	result := map[string]string{}
	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			panic(err)
		}
		result[record[colCode]] = record[colTitle]
	}
	fmt.Println("Loaded names for ", len(result), " ICD codes.")
	return result
}

// LoadICD9ToICD10 loads an ICD-9 to ICD-10 mapping. Supports CSV (.csv) and JSON formats.
func LoadICD9ToICD10(file string) map[string]string {
	ext := strings.ToLower(filepath.Ext(file))
	if ext == ".csv" {
		return loadICD9ToICD10FromCSV(file)
	}
	return loadICD9ToICD10FromJSON(file)
}

func loadICD9ToICD10FromJSON(file string) map[string]string {
	fmt.Println("Loading ICD-9 to ICD-10 mapping from JSON: ", file)
	jsonFile, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer jsonFile.Close()
	jsonBytes, _ := ioutil.ReadAll(jsonFile)
	var mapping map[string]string
	json.Unmarshal(jsonBytes, &mapping)
	fmt.Println("Loaded ", len(mapping), " ICD-9 to ICD-10 mappings.")
	return mapping
}

func loadICD9ToICD10FromCSV(file string) map[string]string {
	fmt.Println("Loading ICD-9 to ICD-10 mapping from CSV: ", file)
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	reader := csv.NewReader(f)
	header, err := reader.Read()
	if err != nil {
		panic(err)
	}
	colICD9, colICD10 := -1, -1
	for i, h := range header {
		switch strings.TrimSpace(strings.ToLower(h)) {
		case "icd9cm":
			colICD9 = i
		case "icd10cm":
			colICD10 = i
		}
	}
	if colICD9 == -1 || colICD10 == -1 {
		panic("ICD9-to-ICD10 CSV must contain icd9cm and icd10cm columns")
	}
	mapping := map[string]string{}
	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			panic(err)
		}
		icd9 := record[colICD9]
		icd10 := record[colICD10]
		if _, exists := mapping[icd9]; !exists {
			mapping[icd9] = icd10
		}
	}
	fmt.Println("Loaded ", len(mapping), " ICD-9 to ICD-10 mappings.")
	return mapping
}

// StripDot removes dots from an ICD code.
func StripDot(code string) string {
	return strings.ReplaceAll(code, ".", "")
}

// ConvertICD9ToICD10 attempts to convert an ICD-9 code using the provided mapping.
// Tries the raw code first, then a dotted variant (dot after 3rd char).
func ConvertICD9ToICD10(icd9Code string, mapping map[string]string) (string, bool) {
	if icd10, ok := mapping[icd9Code]; ok {
		return StripDot(icd10), true
	}
	if len(icd9Code) > 3 {
		dotted := icd9Code[:3] + "." + icd9Code[3:]
		if icd10, ok := mapping[dotted]; ok {
			return StripDot(icd10), true
		}
	}
	return "", false
}

// CheckICD9Coverage checks ICD-9 conversion coverage against diagnoses_icd.csv.
func CheckICD9Coverage(diagnosesICDFile, icd9ToIcd10File string) {
	mapping := LoadICD9ToICD10(icd9ToIcd10File)
	fmt.Println("Scanning ICD-9 codes in: ", diagnosesICDFile)
	file, err := os.Open(diagnosesICDFile)
	if err != nil {
		panic(err)
	}
	defer file.Close()
	reader := csv.NewReader(file)
	header, err := reader.Read()
	if err != nil {
		panic(err)
	}
	colICDCode, colICDVersion := -1, -1
	for i, h := range header {
		switch strings.TrimSpace(strings.ToLower(h)) {
		case "icd_code":
			colICDCode = i
		case "icd_version":
			colICDVersion = i
		}
	}
	if colICDCode == -1 || colICDVersion == -1 {
		panic("diagnoses_icd.csv must contain icd_code and icd_version columns")
	}
	uniqueICD9 := map[string]int{}
	uniqueICD10 := map[string]int{}
	totalRows := 0
	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			panic(err)
		}
		totalRows++
		if record[colICDVersion] == "9" {
			uniqueICD9[record[colICDCode]]++
		} else {
			uniqueICD10[record[colICDCode]]++
		}
	}
	convertible := map[string]bool{}
	unconvertible := map[string]int{}
	for code, count := range uniqueICD9 {
		if _, ok := ConvertICD9ToICD10(code, mapping); ok {
			convertible[code] = true
		} else {
			unconvertible[code] = count
		}
	}
	sumCounts := func(m map[string]int) int {
		total := 0
		for _, c := range m {
			total += c
		}
		return total
	}
	fmt.Println("=== ICD-9 Coverage Report ===")
	fmt.Printf("Total rows in diagnoses_icd.csv: %d\n", totalRows)
	fmt.Printf("  ICD-10 rows: %d (across %d unique codes)\n", sumCounts(uniqueICD10), len(uniqueICD10))
	fmt.Printf("  ICD-9 rows:  %d (across %d unique codes)\n", sumCounts(uniqueICD9), len(uniqueICD9))
	fmt.Printf("ICD-9 codes convertible:     %d / %d\n", len(convertible), len(uniqueICD9))
	fmt.Printf("ICD-9 codes NOT convertible: %d / %d\n", len(unconvertible), len(uniqueICD9))
	lostRows := 0
	for _, count := range unconvertible {
		lostRows += count
	}
	fmt.Printf("ICD-9 rows that would be skipped: %d / %d\n", lostRows, sumCounts(uniqueICD9))
	if len(unconvertible) > 0 {
		fmt.Println("\nUnconvertible ICD-9 codes (up to 30, sorted by frequency):")
		type codeCount struct {
			code  string
			count int
		}
		sorted := make([]codeCount, 0, len(unconvertible))
		for code, count := range unconvertible {
			sorted = append(sorted, codeCount{code, count})
		}
		sort.Slice(sorted, func(i, j int) bool {
			return sorted[i].count > sorted[j].count
		})
		limit := len(sorted)
		if limit > 30 {
			limit = 30
		}
		for i := 0; i < limit; i++ {
			fmt.Printf("  %-10s  (%d occurrences)\n", sorted[i].code, sorted[i].count)
		}
	}
}

// ParseDiagnosesICDWithPatients builds a trajectory.Experiment from diagnoses_icd.csv.
func ParseDiagnosesICDWithPatients(diagnosesICDFile, patientsFile, admissionsFile,
	icdNamesFile, icd9ToIcd10File string, nofCohortAges int,
	filters []trajectory.PatientFilter) (*trajectory.Experiment, *trajectory.PatientMap) {

	demoMap := LoadMIMIC4PatientDemographics(patientsFile)
	admissionDates := LoadAdmissionDates(admissionsFile)
	icdNames := LoadICDNames(icdNamesFile)
	icd9ToIcd10Map := map[string]string{}
	if icd9ToIcd10File != "" {
		icd9ToIcd10Map = LoadICD9ToICD10(icd9ToIcd10File)
	}

	fmt.Println("Parsing MIMIC-IV ICD diagnoses from: ", diagnosesICDFile)
	file, err := os.Open(diagnosesICDFile)
	if err != nil {
		panic(err)
	}
	defer func() {
		if err := file.Close(); err != nil {
			panic(err)
		}
	}()
	reader := csv.NewReader(file)
	header, err := reader.Read()
	if err != nil {
		panic(fmt.Sprint("failed to read diagnoses_icd CSV header: ", err))
	}
	colSubjectID, colHadmID, colICDCode, colICDVersion := -1, -1, -1, -1
	for i, h := range header {
		switch strings.TrimSpace(strings.ToLower(h)) {
		case "subject_id":
			colSubjectID = i
		case "hadm_id":
			colHadmID = i
		case "icd_code":
			colICDCode = i
		case "icd_version":
			colICDVersion = i
		}
	}
	if colSubjectID == -1 || colHadmID == -1 || colICDCode == -1 || colICDVersion == -1 {
		panic("diagnoses_icd.csv must contain subject_id, hadm_id, icd_code, and icd_version columns")
	}

	patientMap := &trajectory.PatientMap{
		PIDMap:       map[int]*trajectory.Patient{},
		PIDStringMap: map[string]int{},
	}
	codeToDID := map[string]int{}
	didToName := map[int]string{}
	didToCode := map[int]string{}
	didCounter := 0
	maxYOB := 1850
	minYOB := 2100
	eventCtr := 0
	skippedNoDate := 0
	skippedNoConvert := 0
	icd9ConvertedCtr := 0
	icd10DirectCtr := 0

	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			panic(err)
		}
		eventCtr++
		subjectIDStr := record[colSubjectID]
		hadmID := record[colHadmID]
		icdCode := record[colICDCode]
		icdVersion := record[colICDVersion]

		date, dateOk := admissionDates[hadmID]
		if !dateOk {
			skippedNoDate++
			continue
		}
		if icdVersion == "9" {
			converted, ok := ConvertICD9ToICD10(icdCode, icd9ToIcd10Map)
			if !ok {
				skippedNoConvert++
				continue
			}
			icd9ConvertedCtr++
			icdCode = converted
		} else {
			icd10DirectCtr++
		}
		icdName := icdCode
		if title, ok := icdNames[icdCode]; ok {
			icdName = title
		}
		if _, exists := codeToDID[icdCode]; !exists {
			codeToDID[icdCode] = didCounter
			didToName[didCounter] = icdName
			didToCode[didCounter] = icdCode
			didCounter++
		}
		if _, exists := patientMap.PIDStringMap[subjectIDStr]; !exists {
			patientMap.Ctr++
			pid := patientMap.Ctr
			sex := trajectory.Male
			yob := 0
			if demo, found := demoMap[subjectIDStr]; found {
				sex = demo.Sex
				yob = demo.Yob
			}
			if sex == trajectory.Male {
				patientMap.MaleCtr++
			} else {
				patientMap.FemaleCtr++
			}
			if yob > 0 {
				maxYOB = utils.MaxInt(yob, maxYOB)
				minYOB = utils.MinInt(yob, minYOB)
			}
			patient := &trajectory.Patient{
				PID: pid, PIDString: subjectIDStr, YOB: yob,
				CohortAge: 0, Sex: sex,
				Diagnoses: []*trajectory.Diagnosis{}, Region: 0,
			}
			patientMap.PIDMap[pid] = patient
			patientMap.PIDStringMap[subjectIDStr] = pid
		}
		pid := patientMap.PIDStringMap[subjectIDStr]
		patient := patientMap.PIDMap[pid]
		did := codeToDID[icdCode]
		diagnosis := &trajectory.Diagnosis{PID: pid, DID: did, Date: date}
		trajectory.AddDiagnosis(patient, diagnosis)
	}

	fmt.Println("Parsed ", eventCtr, " ICD diagnosis rows for ", len(patientMap.PIDMap), " patients.")
	fmt.Printf("Males: %d, Females: %d\n", patientMap.MaleCtr, patientMap.FemaleCtr)
	fmt.Println("ICD-10 direct: ", icd10DirectCtr, ", ICD-9 converted: ", icd9ConvertedCtr)
	fmt.Println("Unique ICD-10 codes: ", didCounter)
	fmt.Println("Skipped (no admission date): ", skippedNoDate)
	fmt.Println("Skipped (ICD-9 not convertible): ", skippedNoConvert)
	fmt.Printf("Year of birth range: %d - %d\n", minYOB, maxYOB)

	return BuildMIMIC4Experiment(patientMap, didToName, didToCode, didCounter,
		nofCohortAges, maxYOB, minYOB, "mimic4_diagnoses_icd", filters)
}

// ---------------------------------------------------------------------------
// Shared experiment builder
// ---------------------------------------------------------------------------

// BuildMIMIC4Experiment computes cohort ages, sorts/compacts diagnoses, applies filters,
// and assembles a trajectory.Experiment from the parsed patient and code data.
func BuildMIMIC4Experiment(patientMap *trajectory.PatientMap,
	didToName map[int]string, didToCode map[int]string, nofDiagnosisCodes,
	nofCohortAges, maxYOB, minYOB int, name string,
	filters []trajectory.PatientFilter) (*trajectory.Experiment, *trajectory.PatientMap) {

	if nofCohortAges > 1 && maxYOB > minYOB {
		ageRange := math.Ceil(float64(maxYOB-minYOB) / float64(nofCohortAges))
		for _, p := range patientMap.PIDMap {
			if p.YOB > 0 {
				p.CohortAge = int(math.Floor(float64(p.YOB-minYOB) / ageRange))
			}
		}
	}
	for _, patient := range patientMap.PIDMap {
		trajectory.SortDiagnoses(patient)
		trajectory.CompactDiagnoses(patient)
	}
	patients := trajectory.ApplyPatientFilters(filters, patientMap)
	fmt.Println("After filtering: ", len(patients.PIDMap), " patients remain.")

	nofRegions := 1
	cohorts := trajectory.InitializeCohorts(patients, nofCohortAges, nofRegions, nofDiagnosisCodes)
	mergedCohort := trajectory.MergeCohorts(cohorts)

	exp := &trajectory.Experiment{
		NofAgeGroups:      nofCohortAges,
		NofRegions:        nofRegions,
		Level:             0,
		NofDiagnosisCodes: nofDiagnosisCodes,
		DxDRR:             trajectory.MakeDxDRR(nofDiagnosisCodes),
		DxDPatients:       trajectory.MakeDxDPatients(nofDiagnosisCodes),
		DPatients:         mergedCohort.DPatients,
		Cohorts:           cohorts,
		Name:              name,
		NameMap:           didToName,
		IdMap:             didToCode,
		MCtr:              patients.MaleCtr,
		FCtr:              patients.FemaleCtr,
	}
	return exp, patients
}
