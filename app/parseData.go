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
	"fmt"
	"io"
	"log/slog"
	"math"
	"os"
	"path/filepath"
	"ptra/trajectory"
	"ptra/utils"
	"strconv"
)

//Generic ptra parsing

// getIcd10DescToExcludeFromAnalysis returns a map that lists ICD10 categories to be excluded from analysis by mapping
// the ICD10 category description (string) onto a boolean.
func getIcd10DescToExcludeFromAnalysis() map[string]bool {
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
func getIcd10CodesToExcludeFromAnalysis() map[string]bool {
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

// Parsing ATC codes
// parseATCCodes parses a files with ATC codes + descriptions. Input: a csv file with the header:
// atc_code,atc_name,ddd,uom,adm_r,note. Returns a map: code -> description.
func parseATCCodes(file string) map[string]string {
	dct := map[string]string{}
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
	reader := csv.NewReader(csvFile)
	reader.Read() //skip header
	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			panic(err)
		}
		dct[record[0]] = record[1]
	}
	return dct
}

// Parsing patient information.
// parsePatientData parses a file with patient information. Input: a patient file in csv
// format, a desired number of age groups to initialize cohorts. Diagnoses of the patient need to be filled in after
// parsing the diagnoses file.
func parsePatientData(file string, nofCohortAges int) (*trajectory.PatientMap, int) {
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
	//the header should be: patient_id, sex, race, ethnicity, year_of_birth,
	//age_at_death, patient_regional_location, postal_code, marital_status, reason_yob_missing, month_year_death,
	//source_id
	//skip header
	reader.Read()
	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			panic(err)
		}
		var yob int
		if yob, err = strconv.Atoi(record[2]); err != nil {
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
		dateOfDeathString := record[4]
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
		region := record[3]
		if _, ok := regions[region]; !ok {
			regions[region] = 0
			regionIds[region] = len(regionIds)
		} else {
			regions[region]++
		}
		control := false
		if record[5] == "yes" {
			control = true
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
			Control:   control,
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

// parseDiagnosisDate turns a date string into DiagnosisDate object.
func parseDiagnosisDate(date string) trajectory.DiagnosisDate {
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

// EventOfInterest checks if the ICD10 code is related to bladder cancer
func EventOfInterest(icd10ID string, options []string) bool {
	for _, option := range options {
		if icd10ID == option { //specific ICD10
			return true
		}
		if len(option) == 3 { //the whole category
			if icd10ID[0:3] == option[0:3] {
				return true
			}
		}
	}
	return false
}

// parsePatientDiagnoses parses a csv file containing patient diagnoses. It fills in those diagnoses for the given
// patients. It uses the icd10AnalysisMap to assign internal analysis DID to the diagnoses.
func parsePatientDiagnoses(diagnosesFile string, patients *trajectory.PatientMap,
	icd10AnalysisMap AnalysisMaps, eoid []string) {
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
	ctrExcl := 0
	EOICtr := 0
	reader.Read() //skip header
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
		DIDString := record[2]
		date := parseDiagnosisDate(record[3])

		nr := icd10AnalysisMap.fillInPatientDiagnoses(patient, DIDString, date)
		if nr > 0 {
			ctrExcl++
			continue
		}
		//Check if diagnosis is event of interest.
		if patient.EOIDate == nil && EventOfInterest(DIDString, eoid) {
			EOICtr++
			patient.EOIDate = &date // mark first event of interest
		}
	}
	for _, patient := range patients.PIDMap {
		trajectory.SortDiagnoses(patient)
		trajectory.CompactDiagnoses(patient)
	}
	slog.Info("Parsed diagnosis data",
		slog.Int("diagnoses", ctr),
		slog.Int("excluded", ctrExcl),
		slog.Int("events-of-interest", EOICtr),
	)
}

func (maps *icd10AnalysisMapsFromXML) extendWithATC(file string) {
	codes := parseATCCodes(file)
	for code, desc := range codes {
		maps.NofDiagnosisCodes++
		maps.DIDMap[code] = maps.NofDiagnosisCodes
		maps.NameMap[maps.NofDiagnosisCodes] = desc
	}
}

func (maps *icd10AnalysisMapsFromCCSR) extendWithATC(file string) {
	codes := parseATCCodes(file)
	for code, desc := range codes {
		maps.NofDiagnosisCodes++
		maps.DIDMap[code] = []int{maps.NofDiagnosisCodes}
		maps.NameMap[maps.NofDiagnosisCodes] = desc
	}
}

func ParseData(name, patientFile, diagnosisFile, diagnosisInfoFile, actFile string, nofCohortAges, level int,
	filters []trajectory.PatientFilter, eoid []string) (*trajectory.Experiment, *trajectory.PatientMap) {
	// parse data
	// fill in patients
	patients, nofRegions := parsePatientData(patientFile, nofCohortAges)
	// fill in icd10 to analysis map
	var analysisMaps AnalysisMaps
	var nofDiagnosisCodes int
	var nameMap map[int]string
	var idMap map[int]string
	if filepath.Ext(diagnosisInfoFile) == ".xml" {
		maps := initializeIcd10AnalysisMapsFromXML(diagnosisInfoFile, level, getIcd10DescToExcludeFromAnalysis())
		if actFile != "" {
			maps.extendWithATC(actFile)
		}
		analysisMaps = maps
		nofDiagnosisCodes = maps.NofDiagnosisCodes
		nameMap = maps.NameMap
		idMap = maps.getIdMap()
	}
	if filepath.Ext(diagnosisInfoFile) == ".csv" || filepath.Ext(diagnosisInfoFile) == ".CSV" {
		maps := initializeIcd10AnalysisMapsFromCCSR(diagnosisInfoFile, getIcd10CodesToExcludeFromAnalysis())
		if actFile != "" {
			maps.extendWithATC(actFile)
		}
		analysisMaps = maps
		nofDiagnosisCodes = maps.NofDiagnosisCodes
		nameMap = maps.NameMap
		idMap = maps.getIdMap()
	}
	// fill in diagnoses for patients
	parsePatientDiagnoses(diagnosisFile, patients, analysisMaps, eoid)
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
		DxDRRPval:         trajectory.MakeDxDRR(nofDiagnosisCodes),
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
