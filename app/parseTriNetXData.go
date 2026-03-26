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
	"bytes"
	"encoding/csv"
	"encoding/json"
	"fmt"
	"io"
	"log/slog"
	"math"
	"os"
	"path/filepath"
	"ptra/trajectory"
	"ptra/utils"
	"sort"
	"strconv"
	"strings"
)

// getIcd10DescToExcludeFromAnalysis returns a map that lists ICD10 categories to be excluded from analysis by mapping
// the ICD10 category description (string) onto a boolean.
func getIcd10DescToExcludeFromTriNetXAnalysis() map[string]bool {
	exclude := map[string]bool{}
	exclude["Pregnancy, childbirth and the puerperium (O00-O9A)"] = true
	exclude["Certain conditions originating in the perinatal period (P00-P96)"] = true
	exclude["Symptoms, signs and abnormal clinical and laboratory findings, not elsewhere classified (R00-R99)"] = true
	exclude["Injury, poisoning and certain other consequences of external causes (S00-T88)"] = true
	exclude["External causes of morbidity (V00-Y99)"] = true
	exclude["Factors influencing health status and contact with health services (Z00-Z99)"] = true
	return exclude
}

// getIcd10CodesToExcludeFromAnalysis returns the first letters of ICD10 codes to exclude from analysis.
func getIcd10CodesToExcludeFromTriNetXAnalysis() map[string]bool {
	exclude := map[string]bool{}
	exclude["O"] = true
	exclude["P"] = true
	exclude["R"] = true
	exclude["S"] = true
	exclude["T"] = true
	exclude["V"] = true
	exclude["X"] = true
	exclude["Y"] = true
	exclude["Z"] = true
	return exclude
}

// getNonICD10CodesToAddToAnalysis returns a set of mockup ICD10 codes to be able to introduce non ICD codes to be
// included for analysis. It returns a map from mockup ICD10 code (string) to description string. It introduces "C98" for
// "Radical custectomy (bladder cancer)", "C99" for "MVAC Chemotherapy (bladder cancer)", and "C100" for "Intravesical
// therapy (bladder cancer)".
func getNonICD10CodesToAddToAnalysis() map[string]string {
	return map[string]string{
		"C98":  "Radical cystectomy (bladder cancer)",
		"C99":  "MVAC Chemotherapy (bladder cancer)",
		"C100": "Intravesical therapy (bladder cancer)",
	}
}

// Parsing patient information.
// parseTriNetXPatientData parses a file with patient information from the TriNetX database. Input: a patient file in csv
// format, a desired number of age groups to initialize cohorts. Diagnoses of the patient need to be filled in after
// parsing the diagnoses file.
func parseTriNetXPatientData(file string, nofCohortAges int) (*trajectory.PatientMap, int) {
	//open file
	csvFile, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer func() {
		if err := csvFile.Close(); err != nil {
			panic(err)
		}
	}()
	patientMap := &trajectory.PatientMap{PIDMap: map[int]*trajectory.Patient{}, PIDStringMap: map[string]int{}}
	maxYOB := 1850
	minYOB := 2021
	deathCr := 0
	regions := map[string]int{} //counts per region
	regionIds := map[string]int{}
	//parse file
	reader := csv.NewReader(csvFile)
	//the header is omitted from the TriNetX file, but is should be: patient_id, sex, race, ethnicity, year_of_birth,
	//age_at_death, patient_regional_location, postal_code, marital_status, reason_yob_missing, month_year_death,
	//source_id
	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			panic(err)
		}
		var yob int
		if yob, err = strconv.Atoi(record[4]); err != nil {
			continue //skip patients without year of birth
		}
		pidString := record[0]
		patientMap.Ctr++      // avoid using 0 as PID
		pid := patientMap.Ctr //analysis ID
		var sex int
		if record[1] == "M" {
			sex = trajectory.Male
			patientMap.MaleCtr++
		}
		if record[1] == "F" {
			sex = trajectory.Female
			patientMap.FemaleCtr++
		}
		dateOfDeathString := record[10]
		var dateOfDeath *trajectory.DiagnosisDate
		if len(dateOfDeathString) == 6 {
			year, err := strconv.Atoi(dateOfDeathString[0:4])
			if err == nil {
				month, err := strconv.Atoi(dateOfDeathString[4:6])
				if err == nil {
					deathCr++
					dateOfDeath = &trajectory.DiagnosisDate{
						Year:  year,
						Month: month,
						Day:   1, //unknown, default to 1
					}
				}
			}
		}
		region := record[6]
		if _, ok := regions[region]; !ok {
			regions[region] = 0
			regionIds[region] = len(regionIds)
		} else {
			regions[region]++
		}
		patient := trajectory.Patient{
			PID:       pid,
			PIDString: pidString,
			YOB:       yob,
			CohortAge: 0,
			Sex:       sex,
			Diagnoses: []*trajectory.Diagnosis{},
			DeathDate: dateOfDeath,
			Region:    regionIds[region],
		}
		patientMap.PIDMap[pid] = &patient
		patientMap.PIDStringMap[pidString] = pid
		maxYOB = utils.MaxInt(yob, maxYOB)
		minYOB = utils.MinInt(yob, minYOB)
	}
	// initialize patient age groups
	ageRange := float64(maxYOB-minYOB) / float64(nofCohortAges)
	ageRange = math.Ceil(ageRange)
	if nofCohortAges > 1 {
		for _, p := range patientMap.PIDMap {
			p.CohortAge = int(math.Floor(float64(p.YOB-minYOB) / float64(ageRange)))
		}
	}

	var buffer bytes.Buffer
	for region, nr := range regions {
		buffer.WriteString(fmt.Sprintf("[%s: %d] ", region, nr))
	}
	slog.Info("Parsed patient data",
		slog.Int("#regions", len(regions)),
	)
	slog.Debug("Parsed patient data",
		slog.String("Regions", buffer.String()),
	)
	slog.Info("Parsed patients with year of birth known",
		slog.Int("N", patientMap.Ctr),
		slog.Int("females", patientMap.FemaleCtr),
		slog.Int("males", patientMap.MaleCtr),
		slog.Int("#deaths", deathCr),
	)
	slog.Info("Year of birth",
		slog.Int("oldest", minYOB),
		slog.Int("youngest", maxYOB),
	)

	return patientMap, len(regions)
}

//Parsing patient diagnoses

// parseTriNetXDiagnosisDate turns a TriNetX date string into DiagnosisDate object.
func parseTriNetXDiagnosisDate(date string) trajectory.DiagnosisDate {
	year, err := strconv.Atoi(date[0:4])
	if err != nil {
		panic(err)
	}
	month, err := strconv.Atoi(date[5:7])
	if err != nil {
		panic(err)
	}
	day, err := strconv.Atoi(date[8:10])
	if err != nil {
		panic(err)
	}
	return trajectory.DiagnosisDate{Year: year, Month: month, Day: day}
}

// TriNetXEventOfInterest checks if the ICD10 code is related to bladder cancer
func TriNetXEventOfInterest(icd10ID string) bool {
	if icd10ID == "Z85.1" {
		return true
	}
	if len(icd10ID) >= 3 && icd10ID[0:3] == "C67" {
		return true
	}
	return false
}

// TreatmentInfo implements a structure for storing the dates of certain bladder cancer treatments.
type TreatmentInfo struct {
	RCDate   *trajectory.DiagnosisDate //Date of radical cystectomy
	MVACDate *trajectory.DiagnosisDate //Date of MVAC chemotherapy
	IVTDate  *trajectory.DiagnosisDate //Date of intravesical therapy
}

// parseTriNetXTreatmentFile parses a csv file that contains information of patient's treatments at different time stamps.
// It returns a map from PID -> TreatmentInfo.
func parseTriNetXTreatmentFile(fileName string) map[string]*TreatmentInfo {
	result := map[string]*TreatmentInfo{}
	file, err := os.Open(fileName)
	if err != nil {
		panic(err)
	}
	defer func() {
		if err := file.Close(); err != nil {
			panic(err)
		}
	}()
	reader := csv.NewReader(file)
	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			panic(err)
		}
		PIDString := record[0]
		var rcDate, mvacDate, ivtDate *trajectory.DiagnosisDate
		if len(record[10]) == 10 { // valid date
			d := parseTriNetXDiagnosisDate(record[10])
			rcDate = &d
		}
		if len(record[11]) == 10 {
			d := parseTriNetXDiagnosisDate(record[11])
			mvacDate = &d
		}
		if len(record[13]) == 10 {
			d := parseTriNetXDiagnosisDate(record[13])
			rcDate = &d
		}
		result[PIDString] = &TreatmentInfo{RCDate: rcDate, MVACDate: mvacDate, IVTDate: ivtDate}
	}
	return result
}

// parseTrinetXPatientDiagnoses parses a csv file containing patient diagnoses. It fills in those diagnoses for the given
// patients. It uses the icd10AnalysisMap to assign internal analysis DID to the diagnoses.
// TO DO: Handle ICD09 diagnoses.
func parseTrinetXPatientDiagnoses(diagnosesFile, treatmentInfoFile string, patients *trajectory.PatientMap, icd10AnalysisMap AnalysisMaps, icd9ToIcd10Map map[string]string) {
	file, err := os.Open(diagnosesFile)
	if err != nil {
		panic(err)
	}
	defer func() {
		if err := file.Close(); err != nil {
			panic(err)
		}
	}()
	reader := csv.NewReader(file)
	ctr := 0 //for counting the number of parsed diagnoses
	ctrID09 := 0
	ctrExcl := 0
	EOICtr := 0
	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			panic(err)
		}
		ctr++
		PIDString := record[0]
		patient, ok := trajectory.GetPatient(PIDString, patients)
		if !ok {
			continue //skip unknown patients
		}
		DIDCodeSystem := record[2]
		DIDString := record[3]
		if (DIDCodeSystem != "ICD-10-CM") && (DIDCodeSystem != "ICD-10-WHO") {
			// try to remap ICD9 code to ICD10 codes
			if DIDString, ok = icd9ToIcd10Map[DIDString]; !ok {
				continue // skip unkown ICD9 codes
			}
			ctrID09++
		}
		date := parseTriNetXDiagnosisDate(record[7])

		nr := icd10AnalysisMap.fillInPatientDiagnoses(patient, DIDString, date)
		if nr > 0 {
			ctrExcl++
			continue
		}
		//Check if diagnosis is event of interest.
		if patient.EOIDate == nil && TriNetXEventOfInterest(DIDString) {
			EOICtr++
			patient.EOIDate = &date // mark first event of interest (e.g. bladder cancers diagnosis)
		}
	}
	var nonICD10DiagnosesMap map[string]*TreatmentInfo
	nonICDCtr := 0
	if treatmentInfoFile != "" {
		nonICD10DiagnosesMap = parseTriNetXTreatmentFile(treatmentInfoFile)
		for _, patient := range patients.PIDMap {
			//fill in non ICD10 diagnoses derived from procedure info
			r := icd10AnalysisMap.fillInNonICDPatientDiagnoses(patient, nonICD10DiagnosesMap)
			nonICDCtr = nonICDCtr + r
		}
	}
	for _, patient := range patients.PIDMap {
		trajectory.SortDiagnoses(patient)
		trajectory.CompactDiagnoses(patient)
	}

	slog.Info("Parsed diagnosis data",
		slog.Int("diagnoses", ctr),
		slog.Int("ICD09", ctrID09),
		slog.Int("ICD10", ctr-ctrID09),
		slog.Int("excluded", ctrExcl),
		slog.Int("events-of-interest", EOICtr),
		slog.Int("non-ICD", nonICDCtr),
	)
}

func ParseTriNetXData(name, patientFile, diagnosisFile, diagnosisInfoFile, treatmentInfoFile string, nofCohortAges,
	level int, icd9ToIcd10File string, filters []trajectory.PatientFilter) (*trajectory.Experiment, *trajectory.PatientMap) {
	// parse data
	// fill in patients
	patients, nofRegions := parseTriNetXPatientData(patientFile, nofCohortAges)
	// fill in icd10 to analysis map
	var analysisMaps AnalysisMaps
	var nofDiagnosisCodes int
	var nameMap map[int]string
	var idMap map[int]string
	if filepath.Ext(diagnosisInfoFile) == ".xml" {
		maps := initializeIcd10AnalysisMapsFromXML(diagnosisInfoFile, level, getIcd10DescToExcludeFromTriNetXAnalysis())
		analysisMaps = maps
		nofDiagnosisCodes = maps.NofDiagnosisCodes
		nameMap = maps.NameMap
		idMap = maps.getIdMap()
	}
	if filepath.Ext(diagnosisInfoFile) == ".csv" || filepath.Ext(diagnosisInfoFile) == ".CSV" {
		maps := initializeIcd10AnalysisMapsFromCCSR(diagnosisInfoFile, getIcd10CodesToExcludeFromTriNetXAnalysis())
		analysisMaps = maps
		nofDiagnosisCodes = maps.NofDiagnosisCodes
		nameMap = maps.NameMap
		idMap = maps.getIdMap()
	}
	icd9ToIcd10Map := map[string]string{}
	if icd9ToIcd10File != "" {
		icd9ToIcd10Map = parseIcd9ToIcd10Mapping(icd9ToIcd10File)
	}
	// fill in diagnoses for patients
	parseTrinetXPatientDiagnoses(diagnosisFile, treatmentInfoFile, patients, analysisMaps, icd9ToIcd10Map)
	// Apply patient filter
	patients = trajectory.ApplyPatientFilters(filters, patients)
	slog.Info("Filtered down to", slog.Int("patients", len(patients.PIDMap)))
	// create cohorts
	cohorts := trajectory.InitializeCohorts(patients, nofCohortAges, nofRegions, nofDiagnosisCodes)
	mergedCohort := trajectory.MergeCohorts(cohorts)
	exp := trajectory.Experiment{
		NofAgeGroups:      nofCohortAges,
		Level:             level,
		NofDiagnosisCodes: nofDiagnosisCodes,
		DxDRR:             trajectory.MakeDxDRR(nofDiagnosisCodes),
		DxDPatients:       trajectory.MakeDxDPatients(nofDiagnosisCodes),
		DPatients:         mergedCohort.DPatients,
		Cohorts:           cohorts,
		Name:              name,
		NameMap:           nameMap,
		NofRegions:        nofRegions,
		IdMap:             idMap,
		FCtr:              patients.FemaleCtr,
		MCtr:              patients.MaleCtr,
	}
	return &exp, patients
}

// opening json file with ICD09 -> ICD10 mapping

func parseIcd9ToIcd10Mapping(file string) map[string]string {
	jsonFile, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer jsonFile.Close()
	slog.Info("Parsing ICD9 to ICD10 mapping from a json file.")
	jsonBytes, _ := io.ReadAll(jsonFile)
	var mapping map[string]string
	json.Unmarshal(jsonBytes, &mapping)
	return mapping
}

// TumorInfo is a struct for storing bladder cancer tumor information concerning: tumor size, tumor lymph nodes, tumor
// metastasis
type TumorInfo struct {
	TStage, NStage, MStage, Stage string
	Date                          trajectory.DiagnosisDate
}

// getTumorStage converts tumor size, number of lymph nodes, and metastatis level into an overall cancer stage.
// T stages: Ta,T1,Tis,T2,T3,T4
// N stages: N0,N1,N2,N3
// M stages: M0,M1
// Stage 0a: Ta,N0,M0
// Stage 0is:Tis,N0,M0 known as carcinoma in situ (CIS)
// Stage I: T1,N0,M0
// Stage II: T2,N0,M0
// Stage IIIA: T3a,T3b, or T4a,N0,M0 --or-- T1 to T4a,N1,M0
// Stage IIIB: T1 to T4a, N2 or N3, M0
// Stage IVA: T4b,any N,M0 or any T, any N, M1a
// Stage IVB: any T, any N, M1b
func getTumorStage(tStage, nStage, mStage string) string {
	if nStage == "N0" && mStage == "M0" {
		switch tStage {
		case "Ta":
			return "0a"
		case "Tis":
			return "0is"
		case "T1":
			return "I"
		case "T2":
			return "II"
		case "T3a", "T3b", "T4a":
			return "IIIA"
		}
	}
	if nStage == "N1" && mStage == "M0" {
		switch tStage {
		case "T1", "T1a", "T1b", "T2", "T2a", "T2b", "T3", "T3a", "T3b", "T4a":
			return "IIIA"
		}
	}
	if (nStage == "N2" || nStage == "N3") && mStage == "M0" {
		switch tStage {
		case "T1", "T1a", "T1b", "T2", "T2a", "T2b", "T3", "T3a", "T3b", "T4", "T4a":
			return "IIIB"
		}
	}
	if tStage == "T4b" && mStage == "M0" {
		return "IVA"
	}
	if mStage == "M1a" {
		return "IVA"
	}
	if mStage == "M1b" {
		return "IVB"
	}
	return tStage + nStage + mStage
}

// tumorIsCISStage checks if tumor is flat or carcinoma in situ (CIS).
func tumorIsCISStage(tumor *TumorInfo) bool {
	return tumor.Stage == "0is"
}

// parsetTriNetXTumorData parses the tumor data from a csv file and returns a map PIDString -> []*TumorInfo.
func ParsetTriNetXTumorData(fileName string) map[string][]*TumorInfo {
	file, err := os.Open(fileName)
	if err != nil {
		panic(err)
	}
	defer func() {
		if err := file.Close(); err != nil {
			panic(err)
		}
	}()
	result := map[string][]*TumorInfo{}
	reader := csv.NewReader(file)
	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			panic(err)
		}
		tumorSite := strings.Split(record[4], ".")
		if tumorSite[0] == "C67" { //only record bladder cancer information
			PIDString := record[0]
			date := parseTriNetXDiagnosisDate(record[1])
			tumorSizeInfo := strings.Split(record[10], "_")
			numberOfLymphNodesInfo := strings.Split(record[11], "_")
			metastaticInfo := strings.Split(record[12], "_")
			if len(tumorSizeInfo) == 1 || len(numberOfLymphNodesInfo) == 1 || len(metastaticInfo) == 1 {
				continue
			}
			tumor := &TumorInfo{Date: date, TStage: tumorSizeInfo[1], NStage: numberOfLymphNodesInfo[1],
				MStage: metastaticInfo[1]}
			tumor.Stage = getTumorStage(tumorSizeInfo[1], numberOfLymphNodesInfo[1], metastaticInfo[1])
			if ts, ok := result[PIDString]; ok {
				result[PIDString] = append(ts, tumor)
			} else {
				result[PIDString] = []*TumorInfo{tumor}
			}
		}
	}
	printTumorInfoSummary(result)
	return result
}

func printTumorInfoSummary(tumorInfo map[string][]*TumorInfo) {
	slog.Info("Parsed tumor info. Found info for", slog.Int("patients", len(tumorInfo)))
	ctr := map[string]int{}
	for _, tumors := range tumorInfo {
		for _, tumor := range tumors {
			ctr[tumor.TStage]++
			ctr[tumor.NStage]++
			ctr[tumor.MStage]++
		}
	}
	stages := []string{}
	for stage, _ := range ctr {
		stages = append(stages, stage)
	}
	sort.Strings(stages)
	for _, stage := range stages {
		slog.Debug("tumor info",
			slog.String("stage", stage),
			slog.Int("entries", ctr[stage]),
		)
	}
}

func printTumorInfo(tumorInfo map[int][]*TumorInfo) {
	slog.Info("Tumor Info", slog.Int("#patients", len(tumorInfo)))
	for pid, infos := range tumorInfo {
		for _, info := range infos {
			slog.Debug("  Patient",
				slog.Int("pid", pid),
				slog.String("TStage", info.TStage),
				slog.String("NStage", info.NStage),
				slog.String("MStage", info.MStage),
				slog.String("Global", info.Stage),
			)
		}
	}
}
