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

// ptra_mimic_iv is a standalone application for patient trajectory analysis on
// MIMIC-IV hospital data.  It reads hcpcsevents.csv (HCPCS procedure events) and
// patients.csv from the MIMIC-IV dataset and runs the full ptra pipeline:
// parse → relative-risk scoring → trajectory building → output → optional clustering.
//
// Usage:
//
//	ptra_mimic_iv <hospDir> <outputPath> [flags]
package main

import (
	"bytes"
	"encoding/csv"
	"encoding/json"
	"flag"
	"fmt"
	"io"
	"io/ioutil"
	"log"
	"math"
	"os"
	"path/filepath"
	"ptra/cluster"
	"ptra/trajectory"
	"ptra/utils"
	"runtime"
	"strconv"
	"strings"
)

const (
	programVersion = 0.1
	programName    = "ptra_mimic_iv"
)

func programMessage() string {
	return fmt.Sprint(programName, " version ", programVersion, " compiled with ", runtime.Version())
}

const helpText = "\nptra_mimic_iv parameters:\n" +
	"ptra_mimic_iv hospDir outputPath\n" +
	"\n" +
	"  hospDir    Path to the MIMIC-IV hosp directory\n" +
	"  outputPath Directory where output files are written\n" +
	"\n" +
	"[--source hcpcs | diagnoses_icd]  Event source (default: hcpcs)\n" +
	"[--ICD9ToICD10File file]          JSON mapping ICD-9 to ICD-10 (for diagnoses_icd)\n" +
	"[--nofAgeGroups nr]\n" +
	"[--maxYears nr]\n" +
	"[--minYears nr]\n" +
	"[--minPatients nr]\n" +
	"[--maxTrajectoryLength nr]\n" +
	"[--minTrajectoryLength nr]\n" +
	"[--name string]\n" +
	"[--iter nr]\n" +
	"[--RR nr]\n" +
	"[--saveRR file]\n" +
	"[--loadRR file]\n" +
	"[--cluster]\n" +
	"[--mclPath string]\n" +
	"[--clusterGranularities list]\n" +
	"[--pfilters male | female | age70+ | age70-]\n" +
	"[--nrOfThreads nr]\n"

// ---------------------------------------------------------------------------
// Flag helpers (same pattern as ptra main.go)
// ---------------------------------------------------------------------------

func parseFlags(flags flag.FlagSet, requiredArgs int, help string) {
	if len(os.Args) < requiredArgs {
		fmt.Fprintln(os.Stderr, "Incorrect number of parameters.")
		fmt.Fprint(os.Stderr, help)
		os.Exit(1)
	}
	flags.SetOutput(ioutil.Discard)
	if err := flags.Parse(os.Args[requiredArgs:]); err != nil {
		x := 0
		if err != flag.ErrHelp {
			fmt.Fprint(os.Stderr, err)
		}
		fmt.Fprint(os.Stderr, help)
		os.Exit(x)
	}
	if flags.NArg() > 0 {
		fmt.Fprint(os.Stderr, "Cannot parse remaining parameters:", flags.Args())
		fmt.Fprint(os.Stderr, help)
		os.Exit(1)
	}
}

func getFileName(s, help string) string {
	switch s {
	case "-h", "--h", "-help", "--help":
		fmt.Fprint(os.Stderr, help)
		os.Exit(1)
	}
	return s
}

// ---------------------------------------------------------------------------
// Patient filter helpers (subset relevant for MIMIC-IV, no tumor info)
// ---------------------------------------------------------------------------

func getPatientFilter(s string) trajectory.PatientFilter {
	id := func(p *trajectory.Patient) bool { return true }
	switch s {
	case "id":
		return id
	case "age70+":
		return trajectory.AboveSeventyAggregator()
	case "age70-":
		return trajectory.LessThanSeventyAggregator()
	case "male":
		return trajectory.FemaleFilter()
	case "female":
		return trajectory.MaleFilter()
	default:
		return id
	}
}

func getPatientFilters(f string) []trajectory.PatientFilter {
	fs := strings.Split(f, ",")
	result := []trajectory.PatientFilter{}
	for _, f := range fs {
		result = append(result, getPatientFilter(f))
	}
	return result
}

// ---------------------------------------------------------------------------
// MIMIC-IV HCPCS parser (streaming)
// ---------------------------------------------------------------------------

// parseMIMIC4Date parses a date string in YYYY-MM-DD format into a trajectory.DiagnosisDate.
func parseMIMIC4Date(dateStr string) (trajectory.DiagnosisDate, bool) {
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

// resolveColumns maps CSV header names to their column indices.
func resolveColumns(header []string) (colSubjectID, colChartDate, colHcpcsCd, colShortDesc int) {
	colSubjectID, colChartDate, colHcpcsCd, colShortDesc = -1, -1, -1, -1
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
	return
}

// parseHCPCSEventsWithPatients parses patients.csv (small, fully loaded) then streams
// hcpcsevents.csv row-by-row to build a trajectory.Experiment.
func parseHCPCSEventsWithPatients(hcpcsFile, patientsFile string, nofCohortAges int,
	filters []trajectory.PatientFilter) (*trajectory.Experiment, *trajectory.PatientMap) {

	// --- load patient demographics into a lookup map ---
	type patientDemo struct {
		sex int
		yob int
	}
	demoMap := map[string]*patientDemo{}

	fmt.Println("Parsing MIMIC-IV patient demographics from: ", patientsFile)
	pFile, err := os.Open(patientsFile)
	if err != nil {
		panic(err)
	}
	pReader := csv.NewReader(pFile)
	// skip header
	if _, err := pReader.Read(); err != nil {
		panic(err)
	}
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
		demoMap[subjectIDStr] = &patientDemo{sex: sex, yob: yob}
	}
	pFile.Close()
	fmt.Println("Loaded demographics for ", len(demoMap), " patients.")

	// --- stream hcpcsevents.csv ---
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
	colSubjectID, colChartDate, colHcpcsCd, colShortDesc := resolveColumns(header)

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

		date, ok := parseMIMIC4Date(chartDateStr)
		if !ok {
			skippedCtr++
			continue
		}

		// register HCPCS code
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

		// register patient
		if _, exists := patientMap.PIDStringMap[subjectIDStr]; !exists {
			patientMap.Ctr++
			pid := patientMap.Ctr

			sex := trajectory.Male
			yob := 0
			if demo, found := demoMap[subjectIDStr]; found {
				sex = demo.sex
				yob = demo.yob
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
				PID:       pid,
				PIDString: subjectIDStr,
				YOB:       yob,
				CohortAge: 0,
				Sex:       sex,
				Diagnoses: []*trajectory.Diagnosis{},
				Region:    0,
			}
			patientMap.PIDMap[pid] = patient
			patientMap.PIDStringMap[subjectIDStr] = pid
		}

		// append event
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
	fmt.Println("Year of birth range: ", minYOB, " — ", maxYOB)

	// compute cohort ages from YOB range
	if nofCohortAges > 1 && maxYOB > minYOB {
		ageRange := math.Ceil(float64(maxYOB-minYOB) / float64(nofCohortAges))
		for _, p := range patientMap.PIDMap {
			if p.YOB > 0 {
				p.CohortAge = int(math.Floor(float64(p.YOB-minYOB) / ageRange))
			}
		}
	}

	// sort and compact
	for _, patient := range patientMap.PIDMap {
		trajectory.SortDiagnoses(patient)
		trajectory.CompactDiagnoses(patient)
	}

	// apply filters
	patients := trajectory.ApplyPatientFilters(filters, patientMap)
	fmt.Println("After filtering: ", len(patients.PIDMap), " patients remain.")

	// build cohorts and experiment
	nofDiagnosisCodes := didCounter
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
		Name:              "mimic4_hcpcs",
		NameMap:           didToName,
		IdMap:             didToCode,
		MCtr:              patients.MaleCtr,
		FCtr:              patients.FemaleCtr,
	}

	return exp, patients
}

// ---------------------------------------------------------------------------
// diagnoses_icd parser — entirely separate from the HCPCS parser
// ---------------------------------------------------------------------------

// loadAdmissionDates loads admissions.csv and returns a map from hadm_id → DiagnosisDate.
// The admittime column (format "YYYY-MM-DD HH:MM:SS") is used as the event date.
func loadAdmissionDates(admissionsFile string) map[string]trajectory.DiagnosisDate {
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
		// admittime format: "2180-05-06 22:23:00" — parse the date portion
		if len(record[colAdmitTime]) >= 10 {
			if date, ok := parseMIMIC4Date(record[colAdmitTime][:10]); ok {
				result[record[colHadmID]] = date
			}
		}
	}
	fmt.Println("Loaded admission dates for ", len(result), " admissions.")
	return result
}

// loadICDNames loads d_icd_diagnoses.csv and returns a map from icd_code → long_title.
func loadICDNames(file string) map[string]string {
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

// loadICD9ToICD10 loads a JSON file mapping ICD-9 codes to ICD-10 codes.
// Format: {"ICD9_CODE": "ICD10_CODE", ...} — same as the TriNetX --ICD9ToICD10File.
func loadICD9ToICD10(file string) map[string]string {
	fmt.Println("Loading ICD-9 to ICD-10 mapping from: ", file)
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

// stripDot removes dots from an ICD code. MIMIC-IV stores codes without dots (e.g. "A001")
// while some mapping files use dotted notation (e.g. "A00.1").
func stripDot(code string) string {
	return strings.ReplaceAll(code, ".", "")
}

// convertICD9ToICD10 attempts to convert an ICD-9 code using the provided mapping.
// It tries the raw code first, then a dotted variant (dot after 3rd char), since mapping
// files may use either format. Returns the ICD-10 code (dot-stripped) and true on success.
func convertICD9ToICD10(icd9Code string, mapping map[string]string) (string, bool) {
	if icd10, ok := mapping[icd9Code]; ok {
		return stripDot(icd10), true
	}
	if len(icd9Code) > 3 {
		dotted := icd9Code[:3] + "." + icd9Code[3:]
		if icd10, ok := mapping[dotted]; ok {
			return stripDot(icd10), true
		}
	}
	return "", false
}

// parseDiagnosesICDWithPatients builds a trajectory.Experiment from diagnoses_icd.csv.
// It loads patient demographics from patients.csv, admission dates from admissions.csv,
// and human-readable names from d_icd_diagnoses.csv. ICD-9 codes are converted to ICD-10
// using the provided mapping; rows that cannot be converted are skipped (same as TriNetX).
//
// The file is streamed row-by-row so that large files can be processed efficiently.
func parseDiagnosesICDWithPatients(diagnosesICDFile, patientsFile, admissionsFile,
	icdNamesFile, icd9ToIcd10File string, nofCohortAges int,
	filters []trajectory.PatientFilter) (*trajectory.Experiment, *trajectory.PatientMap) {

	// --- load patient demographics ---
	type patientDemo struct {
		sex int
		yob int
	}
	demoMap := map[string]*patientDemo{}

	fmt.Println("Parsing MIMIC-IV patient demographics from: ", patientsFile)
	pFile, err := os.Open(patientsFile)
	if err != nil {
		panic(err)
	}
	pReader := csv.NewReader(pFile)
	if _, err := pReader.Read(); err != nil {
		panic(err)
	}
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
		demoMap[subjectIDStr] = &patientDemo{sex: sex, yob: yob}
	}
	pFile.Close()
	fmt.Println("Loaded demographics for ", len(demoMap), " patients.")

	// --- load auxiliary data ---
	admissionDates := loadAdmissionDates(admissionsFile)
	icdNames := loadICDNames(icdNamesFile)

	icd9ToIcd10Map := map[string]string{}
	if icd9ToIcd10File != "" {
		icd9ToIcd10Map = loadICD9ToICD10(icd9ToIcd10File)
	}

	// --- stream diagnoses_icd.csv ---
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

		// look up admission date
		date, dateOk := admissionDates[hadmID]
		if !dateOk {
			skippedNoDate++
			continue
		}

		// handle ICD version: convert ICD-9 to ICD-10
		if icdVersion == "9" {
			converted, ok := convertICD9ToICD10(icdCode, icd9ToIcd10Map)
			if !ok {
				skippedNoConvert++
				continue
			}
			icd9ConvertedCtr++
			icdCode = converted
		} else {
			icd10DirectCtr++
		}

		// look up human-readable name
		icdName := icdCode
		if title, ok := icdNames[icdCode]; ok {
			icdName = title
		}

		// register code
		if _, exists := codeToDID[icdCode]; !exists {
			codeToDID[icdCode] = didCounter
			didToName[didCounter] = icdName
			didToCode[didCounter] = icdCode
			didCounter++
		}

		// register patient
		if _, exists := patientMap.PIDStringMap[subjectIDStr]; !exists {
			patientMap.Ctr++
			pid := patientMap.Ctr

			sex := trajectory.Male
			yob := 0
			if demo, found := demoMap[subjectIDStr]; found {
				sex = demo.sex
				yob = demo.yob
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
				PID:       pid,
				PIDString: subjectIDStr,
				YOB:       yob,
				CohortAge: 0,
				Sex:       sex,
				Diagnoses: []*trajectory.Diagnosis{},
				Region:    0,
			}
			patientMap.PIDMap[pid] = patient
			patientMap.PIDStringMap[subjectIDStr] = pid
		}

		// append diagnosis
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
	fmt.Println("Year of birth range: ", minYOB, " — ", maxYOB)

	// compute cohort ages
	if nofCohortAges > 1 && maxYOB > minYOB {
		ageRange := math.Ceil(float64(maxYOB-minYOB) / float64(nofCohortAges))
		for _, p := range patientMap.PIDMap {
			if p.YOB > 0 {
				p.CohortAge = int(math.Floor(float64(p.YOB-minYOB) / ageRange))
			}
		}
	}

	// sort and compact
	for _, patient := range patientMap.PIDMap {
		trajectory.SortDiagnoses(patient)
		trajectory.CompactDiagnoses(patient)
	}

	// apply filters
	patients := trajectory.ApplyPatientFilters(filters, patientMap)
	fmt.Println("After filtering: ", len(patients.PIDMap), " patients remain.")

	// build cohorts and experiment
	nofDiagnosisCodes := didCounter
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
		Name:              "mimic4_diagnoses_icd",
		NameMap:           didToName,
		IdMap:             didToCode,
		MCtr:              patients.MaleCtr,
		FCtr:              patients.FemaleCtr,
	}

	return exp, patients
}

// ---------------------------------------------------------------------------
// main — full pipeline
// ---------------------------------------------------------------------------

func main() {
	var (
		// required positional arguments
		hospDir    string
		outputPath string
		// optional flags
		nofAgeGroups         int
		maxYears             float64
		minYears             float64
		minPatients          int
		maxTrajectoryLength  int
		minTrajectoryLength  int
		name                 string
		iter                 int
		rr                   float64
		saveRR               string
		loadRR               string
		clust                bool
		mclPath              string
		clusterGranularities string
		pfilters             string
		nrOfThreads          int
		source               string
		icd9ToIcd10File      string
	)

	var flags flag.FlagSet
	flags.IntVar(&nofAgeGroups, "nofAgeGroups", 10, "Number of age groups for cohort stratification.")
	flags.IntVar(&nrOfThreads, "nrOfThreads", 0, "Number of threads (0 = GOMAXPROCS default).")
	flags.Float64Var(&maxYears, "maxYears", 5.0, "Maximum years between consecutive events in a pair.")
	flags.Float64Var(&minYears, "minYears", 0.5, "Minimum years between consecutive events in a pair.")
	flags.IntVar(&minPatients, "minPatients", 50, "Minimum patients for last event in a trajectory.")
	flags.IntVar(&maxTrajectoryLength, "maxTrajectoryLength", 5, "Maximum number of events in a trajectory.")
	flags.IntVar(&minTrajectoryLength, "minTrajectoryLength", 3, "Minimum number of events in a trajectory.")
	flags.StringVar(&name, "name", "mimic4_hcpcs", "Experiment name (used in output filenames).")
	flags.IntVar(&iter, "iter", 400, "Monte Carlo sampling iterations for RR calculation.")
	flags.Float64Var(&rr, "RR", 1.0, "Minimum relative-risk score for considering pairs.")
	flags.StringVar(&saveRR, "saveRR", "", "Save computed RR matrix to this file.")
	flags.StringVar(&loadRR, "loadRR", "", "Load RR matrix from a previous run instead of computing.")
	flags.BoolVar(&clust, "cluster", false, "Cluster trajectories using MCL.")
	flags.StringVar(&mclPath, "mclPath", "", "Path to the mcl binary.")
	flags.StringVar(&clusterGranularities, "clusterGranularities", "40,60,80,100", "MCL granularities (comma-separated).")
	flags.StringVar(&pfilters, "pfilters", "id", "Patient filters: id,male,female,age70+,age70- (comma-separated).")
	flags.StringVar(&source, "source", "hcpcs", "Event source: hcpcs or diagnoses_icd.")
	flags.StringVar(&icd9ToIcd10File, "ICD9ToICD10File", "", "JSON file mapping ICD-9 to ICD-10 codes (for diagnoses_icd source).")

	// parse
	parseFlags(flags, 3, helpText)
	hospDir = getFileName(os.Args[1], helpText)
	outputPath, _ = filepath.Abs(getFileName(os.Args[2], helpText))
	outputPath = outputPath + string(filepath.Separator)

	// derive input file paths from the hosp directory and verify they exist
	patientsFile := filepath.Join(hospDir, "patients.csv")

	missing := false
	if _, err := os.Stat(patientsFile); os.IsNotExist(err) {
		fmt.Fprintf(os.Stderr, "Error: patients.csv not found in %s\n", hospDir)
		missing = true
	}

	// source-dependent file validation
	var hcpcsFile, diagnosesICDFile, admissionsFile, icdNamesFile string
	switch source {
	case "hcpcs":
		hcpcsFile = filepath.Join(hospDir, "hcpcsevents.csv")
		if _, err := os.Stat(hcpcsFile); os.IsNotExist(err) {
			fmt.Fprintf(os.Stderr, "Error: hcpcsevents.csv not found in %s\n", hospDir)
			missing = true
		}
	case "diagnoses_icd":
		diagnosesICDFile = filepath.Join(hospDir, "diagnoses_icd.csv")
		admissionsFile = filepath.Join(hospDir, "admissions.csv")
		icdNamesFile = filepath.Join(hospDir, "d_icd_diagnoses.csv")
		for _, f := range []struct{ path, label string }{
			{diagnosesICDFile, "diagnoses_icd.csv"},
			{admissionsFile, "admissions.csv"},
			{icdNamesFile, "d_icd_diagnoses.csv"},
		} {
			if _, err := os.Stat(f.path); os.IsNotExist(err) {
				fmt.Fprintf(os.Stderr, "Error: %s not found in %s\n", f.label, hospDir)
				missing = true
			}
		}
		if icd9ToIcd10File != "" {
			if _, err := os.Stat(icd9ToIcd10File); os.IsNotExist(err) {
				fmt.Fprintf(os.Stderr, "Error: ICD9ToICD10 mapping file not found: %s\n", icd9ToIcd10File)
				missing = true
			}
		}
	default:
		fmt.Fprintf(os.Stderr, "Error: unknown --source %q (must be hcpcs or diagnoses_icd)\n", source)
		os.Exit(1)
	}

	if missing {
		fmt.Fprintf(os.Stderr, "The hospDir must point to a MIMIC-IV hosp directory containing the required files.\n")
		os.Exit(1)
	}

	// create output directory
	if err := os.MkdirAll(filepath.Dir(outputPath), 0700); err != nil {
		panic(err)
	}

	// log the command
	var command bytes.Buffer
	fmt.Fprint(&command, os.Args[0], " ", hospDir, " ", outputPath)
	fmt.Fprint(&command, " --nofAgeGroups ", nofAgeGroups)
	fmt.Fprint(&command, " --maxYears ", maxYears)
	fmt.Fprint(&command, " --minYears ", minYears)
	fmt.Fprint(&command, " --minPatients ", minPatients)
	fmt.Fprint(&command, " --maxTrajectoryLength ", maxTrajectoryLength)
	fmt.Fprint(&command, " --minTrajectoryLength ", minTrajectoryLength)
	fmt.Fprint(&command, " --name ", name)
	fmt.Fprint(&command, " --iter ", iter)
	fmt.Fprint(&command, " --RR ", rr)
	if saveRR != "" {
		fmt.Fprint(&command, " --saveRR ", saveRR)
	}
	if loadRR != "" {
		fmt.Fprint(&command, " --loadRR ", loadRR)
	}
	if clust {
		fmt.Fprint(&command, " --cluster")
		fmt.Fprint(&command, " --mclPath ", mclPath)
		fmt.Fprint(&command, " --clusterGranularities ", clusterGranularities)
	}
	fmt.Fprint(&command, " --pfilters ", pfilters)
	if nrOfThreads > 0 {
		runtime.GOMAXPROCS(nrOfThreads)
		fmt.Fprint(&command, " --nrOfThreads ", nrOfThreads)
	}

	log.Println(programMessage())
	log.Println("Executing command:\n", command.String())

	// 1. Parse inputs into experiment
	var exp *trajectory.Experiment
	var patients *trajectory.PatientMap
	switch source {
	case "hcpcs":
		exp, patients = parseHCPCSEventsWithPatients(hcpcsFile, patientsFile, nofAgeGroups,
			getPatientFilters(pfilters))
	case "diagnoses_icd":
		exp, patients = parseDiagnosesICDWithPatients(diagnosesICDFile, patientsFile,
			admissionsFile, icdNamesFile, icd9ToIcd10File, nofAgeGroups,
			getPatientFilters(pfilters))
	}
	exp.Name = name

	// 2. Initialise relative risk ratios (or load from file)
	if loadRR != "" {
		trajectory.LoadRRMatrix(exp, loadRR)
		trajectory.LoadDxDPatients(exp, patients, fmt.Sprintf("%s.patients.csv", loadRR))
	} else {
		trajectory.InitializeExperimentRelativeRiskRatios(exp, minYears, maxYears, iter)
	}
	if saveRR != "" {
		trajectory.SaveRRMatrix(exp, saveRR)
		trajectory.SaveDxDPatients(exp, fmt.Sprintf("%s.patients.csv", saveRR))
	}
	// assist the gc
	exp.Cohorts = nil
	exp.DPatients = nil

	// 3. Build trajectories
	trajectory.BuildTrajectories(exp, minPatients, maxTrajectoryLength, minTrajectoryLength,
		minYears, maxYears, rr, []trajectory.TrajectoryFilter{})

	// 4. Output trajectories
	trajectory.PrintTrajectoriesToFile(exp, outputPath)
	fmt.Println("Collected trajectories: ")
	for i := 0; i < utils.MinInt(len(exp.Trajectories), 100); i++ {
		trajectory.PrintTrajectory(exp.Trajectories[i], exp)
	}

	// 5. Optional clustering
	if clust {
		var clusterGranularityList []int
		for _, g := range strings.Split(clusterGranularities, ",") {
			gi, _ := strconv.ParseInt(g, 10, 0)
			clusterGranularityList = append(clusterGranularityList, int(gi))
		}
		fmt.Println("MCL Clustering:")
		cluster.ClusterTrajectoriesDirectly(exp, clusterGranularityList, outputPath, mclPath)
	}
}
