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

package lib

import (
	"encoding/csv"
	"fmt"
	"github.com/exascience/pargo/parallel"
	"github.com/imec-int/ptra/lib/utils"
	"github.com/valyala/fastrand"
	"io"
	"math"
	"os"
	"sort"
	"strconv"
	"strings"
)

// Trajectory holds all data relevant to a disease trajectory.
type Trajectory struct {
	Diagnoses      []int            // A list of diagnosis codes that represent the trajectory
	PatientNumbers []int            // A list with nr of patients for each transition in the trajectory
	Patients       [][]*Patient     // A list of patients with the given trajectory
	TrajMap        map[*Patient]int // Maps patient IDs onto a diagnosis index for trajectory tracking
	ID             int              // An analysis id
	Cluster        int              // A cluster ID to which this trajectory is assigned to
}

const (
	Male   = 0
	Female = 1
)

// Patient represents patient information.
type Patient struct {
	PID       int            // analysis ID
	PIDString string         // ID from TriNetX
	YOB       int            // year of birth
	CohortAge int            // age range a patient belongs to
	Sex       int            // 0 = male, 1 = female
	Diagnoses []*Diagnosis   // list of patient's diagnoses, sorted by date <, unique diagnosis per date
	EOIDate   *DiagnosisDate // Event of interest date, e.g. day of cancer diagnosis
	DeathDate *DiagnosisDate // Date of death
	Region    int            // Region where the patient lives
}

// AppendPatient appends a patient to a slice of patients, unless that patient is already a member of that slice.
func AppendPatient(plist []*Patient, p *Patient) []*Patient {
	if !memberPatient(p, plist) {
		return append(plist, p)
	}
	return plist
}

// memberPatient checks if a patient occurs as an entry in a list of patients.
func memberPatient(p *Patient, plist []*Patient) bool {
	for _, p2 := range plist {
		if p2.PID == p.PID {
			return true
		}
	}
	return false
}

// intersectPatients returns the intersection of two lists of patients.
func intersectPatients(patients1 []*Patient, patients2 []*Patient) []*Patient {
	patients3 := []*Patient{}
	for _, p := range patients2 {
		if memberPatient(p, patients1) {
			patients3 = append(patients3, p)
		}
	}
	return patients3
}

// DiagnosisDate represents the date of a diagnosis, with fields for representing the year, month, and day of diagnosis.
type DiagnosisDate struct {
	Year, Month, Day int
}

func DiagnosisDateSmallerThan(d1, d2 DiagnosisDate) bool {
	if d1.Year < d2.Year {
		return true
	}
	if d1.Year > d2.Year {
		return false
	}
	if d1.Month < d2.Month {
		return true
	}
	if d1.Month > d2.Month {
		return false
	}
	if d1.Day < d2.Day {
		return true
	}
	return false
}

// DiagnosisDateToFloat converts a diagnosis date to a floating point number.
func DiagnosisDateToFloat(d DiagnosisDate) float64 {
	return float64(d.Year) + float64(d.Month)/12.0 + float64(d.Day)/365.0
}

// Diagnosis represents a diagnosis for a patient.
type Diagnosis struct {
	PID, DID int
	Date     DiagnosisDate
	Icd10    Icd10Entry
}

// AddDiagnosis appends a diagnosis to a patient's list of diagnoses.
func (p *Patient) AddDiagnosis(d *Diagnosis) {
	p.Diagnoses = append(p.Diagnoses, d)
}

// SortDiagnoses modifies a given patient's list of diagnoses to be ordered by date.
func SortDiagnoses(p *Patient) {
	diagnoses := p.Diagnoses
	sort.Slice(diagnoses, func(i, j int) bool {
		return DiagnosisDateSmallerThan(diagnoses[i].Date, diagnoses[j].Date)
	})
}

// diagnosisDateEqual compares two diagnosis dates for equality in terms of year, month, and day of occurrence.
func diagnosisDateEqual(d1, d2 DiagnosisDate) bool {
	return d1.Year == d2.Year && d1.Month == d2.Month && d1.Day == d2.Day
}

// diagnosisEqual checks if two diagnoses are the same.
func diagnosisEqual(d1, d2 *Diagnosis) bool {
	return d1.DID == d2.DID && d1.PID == d2.PID && diagnosisDateEqual(d1.Date, d2.Date)
}

// CompactDiagnoses makes a sorted diagnosis list contain unique diagnoses for a patient. Want to avoid over counting diagnoses.
func CompactDiagnoses(p *Patient) {
	if len(p.Diagnoses) > 1 {
		diagnoses := p.Diagnoses
		curDiagnosis := diagnoses[0]
		newDiagnoses := []*Diagnosis{}
		newDiagnoses = append(newDiagnoses, curDiagnosis)
		for _, diagnosis := range diagnoses[1:] {
			if !diagnosisEqual(curDiagnosis, diagnosis) {
				curDiagnosis = diagnosis
				newDiagnoses = append(newDiagnoses, curDiagnosis)
			}
		}
		p.Diagnoses = newDiagnoses
	}
}

// PatientMap contains all patient information parsed from the input.
type PatientMap struct {
	PIDStringMap map[string]int   // maps patient string id onto an int PID
	Ctr          int              // total nr of patients parsed, also used for creating PIDs
	PIDMap       map[int]*Patient // maps PID onto a patient object that contain YOB, sex, age group, etc
	// optional info for logging
	MaleCtr   int
	FemaleCtr int
}

// GetPatient retrieves from a patient map the patient object associated with a given patient ID. The patient ID is
// passed as a string and refers to the PID that occurs in the input.
func GetPatient(pidString string, patients *PatientMap) (*Patient, bool) {
	pid, ok := patients.PIDStringMap[pidString]
	if !ok {
		return &Patient{}, false
	}
	patient, ok := patients.PIDMap[pid]
	return patient, ok
}

//Experiment representation

// Cohort represents a specific group of patients from the population stratified by age, sex, and region. The population
// is divided into male and female cohorts. Those cohorts are in turn split into cohorts depending on an age range, e.g.
// this could be one for each possible age range apart by 10 years: [0-10], [10-20],[20-30]...[100-120].
type Cohort struct {
	AgeGroup, Sex, Region, NofPatients, NofDiagnoses int
	DCtr                                             []int        //counts nr of patients per DID
	DPatients                                        [][]*Patient //contains a list of patients per DID
	Patients                                         []*Patient   //the patients in this cohort
}

// MakeDxDRR makes a diagnosis by diagnosis-sized matrix for storing the relative risk score for each possible diagnosis
// pair.
func MakeDxDRR(size int) [][]float64 {
	DxDRR := make([][]float64, size)
	for i := range DxDRR {
		vec := make([]float64, size)
		for j := range vec {
			vec[j] = 1.0
		}
		DxDRR[i] = vec
	}
	return DxDRR
}

// MakeDxDPatients makes a diagnosis by diagnosis-sized matrix for storing the list of patients for each possible
// diagnosis pair.
func MakeDxDPatients(size int) [][][]*Patient {
	DxDPatients := make([][][]*Patient, size)
	for idx := range DxDPatients {
		DxDPatients[idx] = make([][]*Patient, size)
	}
	return DxDPatients
}

// Experiment contains the inputs and outputs for calculating diagnosis trajectories for a specific patient population.
type Experiment struct {
	NofAgeGroups, NofRegions, Level, NofDiagnosisCodes int
	DxDRR                                              [][]float64        // per disease pair, relative risk score (RR)
	DxDPatients                                        [][][]*Patient     // per disease pair, all patients diagnosed
	DPatients                                          [][]*Patient       // per disease, all patients diagnosed
	Cohorts                                            []*Cohort          // cohorts in the experiment
	Name                                               string             // Name of the experiment, for printing
	Icd10Map                                           map[int]Icd10Entry // maps diagnosis ID to Icd10Entry
	Trajectories                                       []*Trajectory      // a list of computed trajectories
	Pairs                                              []*Pair            // a list of all selected pairs that are used to compute trajectories
	IdMap                                              map[int]string     // maps the analysis DID to the original diagnostic ID used in the input data
	MCtr, FCtr                                         int                // counters for counting nr of males,females,patients
}

// selectCohort returns from a list of cohorts a cohort that matches a specific age group, sex, and region.
func selectCohort(cohorts []*Cohort, nofAgeGroups, nofRegions, sex, ageGroup, region int) *Cohort {
	cIndex := cohortIndex(nofAgeGroups, nofRegions, sex, ageGroup, region)
	return cohorts[cIndex]
}

// cohortIndex computes the index of a specific cohort in a cohort array. This index is derived from the sex and age
// group:
// cohorts: [Males: [age: 10-20] [age: 20-30] ... [age: 100-120] Females: [age: 10-20], [age: 20-30] ... [age: 100-120]]
func cohortIndex(nofAgegroups, nofRegions, sex, ageGroup, region int) int {
	return sex*nofAgegroups + ageGroup
}

// makeCohorts creates cohorts for a requested nr of age groups, nr of regions, and nr of diagnosis codes used in patient
// records. Creates empty cohorts for both male and females, for every age group, one for each possible age range.
func makeCohorts(nofAgeGroups, nofRegions, nofDiagnoses int) []*Cohort {
	// Create empty cohorts
	nofCohorts := nofAgeGroups * 2 //#age groups x #sexes
	cohorts := make([]*Cohort, nofCohorts)
	sex := Male
	ageGroup := 0
	region := 0
	femaleCohortIndex := int(math.Floor(float64(nofCohorts / 2)))
	for i := 0; i < nofCohorts; i++ {
		// first fill in male cohorts, then female cohorts
		// check if need to switch to filling in female cohorts
		if i >= femaleCohortIndex && sex != Female {
			sex = Female
			ageGroup = 0
			region = 0
		}
		cohort := &Cohort{
			AgeGroup:     ageGroup,
			Sex:          sex,
			NofPatients:  0,
			NofDiagnoses: 0,
			Region:       region,
			DCtr:         make([]int, nofDiagnoses),
			DPatients:    make([][]*Patient, nofDiagnoses),
			Patients:     []*Patient{},
		}
		cohorts[i] = cohort
		if ageGroup == nofAgeGroups-1 { //switch to next region
			ageGroup = 0
			region++
		} else {
			ageGroup++ //next age group
		}
	}
	return cohorts
}

// InitCohorts creates cohorts + initializes them with the counts for each diagnosis + patients per diagnosis
func InitCohorts(patients *PatientMap, nofAgegroups, nofRegions, nofDiagnosisCodes int) []*Cohort {
	fmt.Println("Initializing cohorts: with ", len(patients.PIDMap), " patients (Males: ", patients.MaleCtr, ""+
		"Females: ", patients.FemaleCtr, ") "+
		" nr of diagnosis codes: ", nofDiagnosisCodes, "nr of age groups: ", nofAgegroups)
	fmt.Println("Making cohort vectors...")
	cohorts := makeCohorts(nofAgegroups, nofRegions, nofDiagnosisCodes)
	// count occurence of diagnoses, collect patients in the cohort
	fmt.Println("Counting diagnosis occurrences...")
	for _, patient := range patients.PIDMap {
		diagnoses := patient.Diagnoses
		cohort := selectCohort(cohorts, nofAgegroups, nofRegions, patient.Sex, patient.CohortAge, patient.Region)
		cohort.NofPatients++
		cohort.Patients = append(cohort.Patients, patient)
		diagnosisCountedForPatient := map[int]bool{} // can count exposure of a disease only once per patient DID->bool
		for _, d1 := range diagnoses {
			// count diagnosis unless already counted (one exposure per patient)
			if _, ok := diagnosisCountedForPatient[d1.DID]; !ok {
				cohort.DCtr[d1.DID]++
				cohort.NofDiagnoses = cohort.NofDiagnoses + 1
				cohort.DPatients[d1.DID] = append(cohort.DPatients[d1.DID], patient)
				diagnosisCountedForPatient[d1.DID] = true
			}
		}
	}
	return cohorts
}

// selectRandomPatientsWithoutShuffle randomly selects number of patients (ctr) from a given list of patients (patients),
// while avoiding patients from a list to be excluded from selection (patientsToExclude). It performs this random selection
// without shuffling the input patients, which would be computationally too costly.
func selectRandomPatientsWithoutShuffle(patients []*Patient, ctr int, patientsToExclude map[int]bool) []*Patient {
	collectedPatients := []*Patient{}
	maxRandSkips := utils.MaxInt(0, len(patients)-len(patientsToExclude)-ctr)
	for _, p := range patients {
		if len(collectedPatients) == ctr {
			break
		}
		if _, ok := patientsToExclude[p.PID]; !ok { // not a member of patients to exclude
			if maxRandSkips > 0 {
				if fastrand.Uint32n(2) > 0 {
					collectedPatients = append(collectedPatients, p)
				} else {
					maxRandSkips--
				}
			} else {
				collectedPatients = append(collectedPatients, p)
			}
		}
	}
	return collectedPatients
}

// selectRandomPatientsFromSimilarCohorts collects for a given list of patients a random list of patients that is
// comparable in terms of cohorts. This means, for each patient, randomly select another patient that belongs to the same
// sex and age groups.
func selectRandomPatientsFromSimilarCohorts(exp *Experiment, patients []*Patient, pids map[int]bool) []*Patient {
	// for each cohort, see how many patients you need to select from it
	cohortSimilar := make([][]*Patient, len(exp.Cohorts))
	for i := range cohortSimilar {
		cohortSimilar[i] = []*Patient{}
	}
	for _, p := range patients {
		cohortIndex := cohortIndex(exp.NofAgeGroups, exp.NofRegions, p.Sex, p.CohortAge, p.Region)
		cohortSimilar[cohortIndex] = append(cohortSimilar[cohortIndex], p)
	}
	// select Random patients from the cohorts
	collectedPatients := []*Patient{}
	for i, ps := range cohortSimilar {
		similarPatients := selectRandomPatientsWithoutShuffle(exp.Cohorts[i].Patients, len(ps), pids)
		for _, p := range similarPatients {
			collectedPatients = append(collectedPatients, p)
		}
	}
	return collectedPatients
}

// probNotExposed calculates for a list of patients exposed to a disease d1, the chance to select a patient exposed to d2
// that is not exposed to d1.
func probNotExposed(exp *Experiment, d1Patients []*Patient, d1IDs map[int]bool, d2 int) float64 {
	d2Ctr := 0.0
	for _, p := range d1Patients {
		idx := cohortIndex(exp.NofAgeGroups, exp.NofRegions, p.Sex, p.CohortAge, p.Region)
		cohort := exp.Cohorts[idx]
		d2Patients := cohort.DPatients[d2]
		ctr := 0
		for _, p2 := range d2Patients { //d2 patients without d1 that could potentially be sampled from
			if _, ok := d1IDs[p2.PID]; !ok { // not a d1 patient
				ctr++
			}
		}
		d2Ctr = d2Ctr + (float64(ctr) / float64(cohort.NofPatients))
	}
	return d2Ctr / float64(len(d1Patients))
}

// countPatientDiagnosis returns 1 if a patient has been diagnosed with a disease (did) or 0 when not.
func countPatientDiagnosis(p *Patient, did int) int {
	for _, d := range p.Diagnoses {
		if d.DID == did {
			return 1
		}
	}
	return 0
}

// countPatientDiagnosisPair returns 1 when a patient was diagnosed with a specific diagnosis pair (d1->d2) and 0 when
// not diagnosed.
func countPatientDiagnosisPair(p *Patient, d1, d2 int, minTime, maxTime float64) (int, int) {
	var d1Date DiagnosisDate
	var d1Index int
	d1ok := false
	for i, d := range p.Diagnoses {
		if d.DID == d1 {
			d1Date = d.Date
			d1Index = i
			d1ok = true
			break
		}
	}
	if d1ok != true {
		panic(fmt.Sprint("Disease d1: ", d1, " not present in patient when checking for d1->d2"))
	}
	for i, d := range p.Diagnoses[d1Index+1:] {
		if d.DID == d2 {
			timeBetween := DiagnosisDateToFloat(d.Date) - DiagnosisDateToFloat(d1Date)
			if timeBetween <= maxTime && timeBetween >= minTime {
				return 1, i
			}
		}
	}
	return 0, -1
}

// countPatientTrajectory returns an index in a patient's diagnosis list when the patient was diagnosed with a diagnosis
// (d) with ond this diagnosis occurs within a specific time frame (cf. minTime and maxTime) of a previous diagnosis
// occuring at index idx in the patient's diagnosis list.
func countPatientTrajectory(p *Patient, idx, d2 int, minTime, maxTime float64) int {
	d1Date := p.Diagnoses[idx].Date
	for i := idx; i < len(p.Diagnoses); i++ {
		diag := p.Diagnoses[i]
		if diag.DID == d2 {
			timeBetween := DiagnosisDateToFloat(diag.Date) - DiagnosisDateToFloat(d1Date)
			if timeBetween <= maxTime && timeBetween >= minTime {
				return i
			}
		}
	}
	return -1
}

// patientsToIdMap creates a map from PID->bool so that for each patient in a given list (patients), there is an entry
// in that map.
func patientsToIdMap(patients []*Patient) map[int]bool {
	pids := map[int]bool{}
	for _, p := range patients {
		pids[p.PID] = true
	}
	return pids
}

// InitRR computes the relative risk ratios for each possible diagnosis pair in an
// experiment. It takes into account the minimum and maximum time between diagnoses (minTime and maxTime). It is an
// iterative algorithm that runs for a given number of iterations (iter). With iter = 400, the calculated p-values are
// within 0.05 of the true p-values and with iter = 10000 they are within 0.01 of the true p-values.
// The relative risk ratios are calculated in parallel for all possible diagnosis pairs.
func (exp *Experiment) InitRR(minTime, maxTime float64, iter int) {
	fmt.Println("Initializing relative risk ratios...")
	fmt.Println("Sampling ", iter, " comparison groups for each diagnosis pair...")
	var indexVector []int
	for i := 0; i < exp.NofDiagnosisCodes; i++ {
		indexVector = append(indexVector, i)
	}
	parallel.Range(0, len(indexVector), 0, func(low, high int) {
		for _, d1 := range indexVector[low:high] {
			d1ExposedPatients := exp.DPatients[d1]
			d1ExposedPatientsIDMap := patientsToIdMap(d1ExposedPatients)
			if len(d1ExposedPatients) > 0 {
				parallel.Range(0, len(indexVector), 0, func(low, high int) {
					for _, d2 := range indexVector[low:high] {
						// select randomly patients without d1 as a control group of same size as group 1
						notd1ExposedPatients := selectRandomPatientsFromSimilarCohorts(exp, d1ExposedPatients, d1ExposedPatientsIDMap)
						if len(d1ExposedPatients) == len(notd1ExposedPatients) {
							// count nr of patients with d2 in the exposed group, taking into account time constraints
							// between exposure and diagnosis d1
							d2CtrInExposedGroup := 0
							d1FollowedByd2Patients := []*Patient{}
							for _, p := range d1ExposedPatients {
								ctr, _ := countPatientDiagnosisPair(p, d1, d2, minTime, maxTime)
								if ctr > 0 {
									d1FollowedByd2Patients = AppendPatient(d1FollowedByd2Patients, p)
								}
								d2CtrInExposedGroup = d2CtrInExposedGroup + ctr
							}
							// count nr of patients with d2 in the not exposed group
							// take the average of this of 400 iterations; 400 iterations to get within 0.05 of the
							// true p-value.
							// first filter out pairs (d1, d2) with a high chance that #d2 in non exposed >= #d1->d2 in exposed
							probd2Notd1Exposed := probNotExposed(exp, d1ExposedPatients, d1ExposedPatientsIDMap, d2)
							probd2d1Exposed := float64(d2CtrInExposedGroup) / float64(len(d1ExposedPatients))
							if probd2Notd1Exposed >= probd2d1Exposed {
								continue // skip sampling for testing d1->d2 pair because it is unlikely
							}
							var pval float64
							d2CtrInNotExposedGroup := 0 // will be average if N iterations
							for i := 0; i < iter; i++ {
								d2Ctr := 0
								for _, p := range notd1ExposedPatients {
									ctr := countPatientDiagnosis(p, d2)
									d2Ctr = d2Ctr + ctr
									d2CtrInNotExposedGroup = d2CtrInNotExposedGroup + ctr
								}
								if d2Ctr >= d2CtrInExposedGroup { // if #D2 in comparison group >= #D1->D2 in exposed group, unlikely that D1->D2
									pval++
								}
								notd1ExposedPatients = selectRandomPatientsFromSimilarCohorts(exp, d1ExposedPatients, d1ExposedPatientsIDMap)
							}
							pval = pval / float64(iter)
							d2CtrInNotExposedGroup = d2CtrInNotExposedGroup / iter // take the average of d2s counted in all sampled non exposed groups
							if pval > 0.001 {
								continue // seems that #D2 in non-exposed > #D1->D2 in exposed, so unlikely D1->D2
							}
							// compute RR
							a := float64(d2CtrInExposedGroup)
							b := float64(len(d1ExposedPatients) - d2CtrInExposedGroup)
							c := float64(d2CtrInNotExposedGroup)
							d := float64(len(d1ExposedPatients) - d2CtrInNotExposedGroup) //take len(d1ExposedPatients) cause we want same length randomly selected groups
							p1 := a / (a + b)
							p2 := c / (c + d)
							RR := p1 / p2
							// initialize RR, d1->d2 ctrs etc
							exp.DxDRR[d1][d2] = RR
							exp.DxDPatients[d1][d2] = d1FollowedByd2Patients
						}
					}
				})
			}
		}
	})
}

// LoadRRMatrix loads an RR matrix from file and stores it in the given experiment. This file was created from a
// previous run. This can be used instead of initializeRelativeRiskRatiosParallel
func (exp *Experiment) LoadRRMatrix(path string) {
	// map icd10 names to DIDs
	nameMap := map[string]int{}
	for did, icd10 := range exp.Icd10Map {
		nameMap[icd10.Name] = did
	}
	file, err := os.Open(path)
	if err != nil {
		panic(err)
	}
	defer func() {
		if err := file.Close(); err != nil {
			panic(err)
		}
	}()
	reader := csv.NewReader(file)
	reader.Comma = '\t'
	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			panic(err)
		}
		d1 := nameMap[record[0]]
		d2 := nameMap[record[1]]
		RR, err := strconv.ParseFloat(record[2], 64)
		if err != nil {
			panic(err)
		}
		exp.DxDRR[d1][d2] = RR
	}
}

// LoadDxDPatients loads the DxD patients from a file created during a previous run. It takes as parameters the experiment
// (exp), the patient map where to look up patient objects parsed from the input (pMap), and the file Name (path). The
// function side effects the experiment's DxDPatients slice. It uses the pMap to match concrete patient objects with IDs
// stored in the file to be able to fill the patients for each pair of diagnoses.
func (exp *Experiment) LoadDxDPatients(pMap *PatientMap, path string) {
	// map icd10 names to DIDs
	nameMap := map[string]int{}
	for did, icd10 := range exp.Icd10Map {
		nameMap[icd10.Name] = did
	}
	//init DxDPatients map
	for i, js := range exp.DxDPatients {
		for j := range js {
			exp.DxDPatients[i][j] = []*Patient{}
		}
	}
	file, err := os.Open(path)
	if err != nil {
		panic(err)
	}
	defer func() {
		if err := file.Close(); err != nil {
			panic(err)
		}
	}()
	reader := csv.NewReader(file)
	reader.Comma = '\t'
	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			panic(err)
		}
		d1 := nameMap[record[0]]
		d2 := nameMap[record[1]]
		pidStrings := strings.Split(record[2], ",")
		for _, pidString := range pidStrings {
			pid := pMap.PIDStringMap[pidString]
			p := pMap.PIDMap[pid]
			exp.DxDPatients[d1][d2] = append(exp.DxDPatients[d1][d2], p)
		}
	}
}

// SaveRRMatrix stores the RR matrix calculated for the given experiment. The diagnosis pairs from the matrix are
// stored line per line as follows: medical Name 1, medical Name 2, RR.
func (exp *Experiment) SaveRRMatrix(path string) {
	file, err := os.Create(path)
	if err != nil {
		panic(err)
	}
	defer func() {
		if err := file.Close(); err != nil {
			panic(err)
		}
	}()
	for i, js := range exp.DxDRR {
		for j, RR := range js {
			fmt.Fprintf(file, "%s\t%s\t%s\n", exp.Icd10Map[i].Name, exp.Icd10Map[j].Name,
				strconv.FormatFloat(RR, 'E', -1, 64))
		}
	}
}

// SaveDxDPatients saves per disease the PIDs that are diagnosed with this disease
func (exp *Experiment) SaveDxDPatients(path string) {
	file, err := os.Create(path)
	if err != nil {
		panic(err)
	}
	defer func() {
		if err := file.Close(); err != nil {
			panic(err)
		}
	}()
	for did, patientPairs := range exp.DxDPatients {
		for pairIdx, pair := range patientPairs {
			if len(patientPairs) > 0 {
				patients := ""
				for patientIdx, patient := range pair {
					patients = patients + patient.PIDString
					if patientIdx != len(pair)-1 {
						patients = patients + ","
					}
				}
				fmt.Fprintf(file, "%s\t%s\t%s\n", exp.Icd10Map[did].Name, exp.Icd10Map[pairIdx].Name, patients)
			}
		}
	}
}

// MergeCohorts returns a single cohort that merges a list of input cohorts. The goal is to obtain a merged list of
// patients, a merged patient total, and a merged disease total.
func MergeCohorts(cohorts []*Cohort) *Cohort {
	cohort1 := cohorts[0]
	for _, cohort2 := range cohorts[1:] {
		// merge patient ctr
		cohort1.NofPatients = cohort1.NofPatients + cohort2.NofPatients
		cohort1.NofDiagnoses = cohort1.NofDiagnoses + cohort2.NofDiagnoses
		// merge DCtr
		for i, ctr := range cohort2.DCtr {
			cohort1.DCtr[i] = cohort1.DCtr[i] + ctr
		}
		// merge DPatients
		for i, ps := range cohort2.DPatients {
			for _, p := range ps {
				cohort1.DPatients[i] = append(cohort1.DPatients[i], p)
			}
		}
	}
	fmt.Println("Merged cohort")
	cohort1.Log(22)
	return cohort1
}

// Log prints a cohort to standard output.
func (cohort Cohort) Log(max int) {
	fmt.Println("Cohort: ")
	fmt.Println("Age group: ", cohort.AgeGroup, " Sex: ", cohort.Sex, " Region: ", cohort.Region, " Nr of patients: ", cohort.NofPatients, " "+
		"Nr of diagnoses: ", cohort.NofDiagnoses)
	fmt.Println("DCtr: [")
	for i := 0; i < utils.MinInt(max, len(cohort.DCtr)); i++ {
		fmt.Print(cohort.DCtr[i], ", ")
	}
	fmt.Println("...]")
}

// Pair is a struct for representing a diagnosis pair. It simply stores two diagnosis codes.
type Pair struct {
	First, Second int
}

// selectDiagnosisPairs selects diagnosis pairs from which to calculate trajectories. These pairs are constrained by
// requiring a minimum number of patients that is diagnosed with the disease pair, and a minimum RR score.
func (exp *Experiment) selectDiagnosisPairs(minPatients int, minRR float64) []*Pair {
	fmt.Println("Selecting diagnosis pairs for building trajectories...")
	var pairs []*Pair
	// fixme: should we use exp.NofDiagnosisCodes?
	nofDiagnosisCodes := len(exp.Icd10Map)
	for i := 0; i < nofDiagnosisCodes; i++ {
		for j := i; j < nofDiagnosisCodes; j++ {
			occurs := len(exp.DxDPatients[i][j])
			occursReverse := len(exp.DxDPatients[j][i])
			RR := exp.DxDRR[i][j]
			RRReverse := exp.DxDRR[j][i]
			if i != j {
				if occurs >= minPatients && RR > minRR && occursReverse >= minPatients && RRReverse > minRR {
					var maxOccurs int
					var maxIndices *Pair
					if occurs > occursReverse {
						maxOccurs = occurs
						maxIndices = &Pair{First: i, Second: j}
					} else {
						maxOccurs = occursReverse
						maxIndices = &Pair{First: j, Second: i}
					}
					test := utils.BinomialCdf(0.5, occurs+occursReverse, maxOccurs)
					if test < 0.05 {
						pairs = append(pairs, maxIndices)
					}
					continue
				}
				if occurs >= minPatients && RR > minRR {
					pairs = append(pairs, &Pair{First: i, Second: j})
					continue
				}
				if occursReverse >= minPatients && RRReverse > minRR {
					pairs = append(pairs, &Pair{First: j, Second: i})
				}
			}
		}
	}
	fmt.Println("Found ", len(pairs), " suitable diagnosis pairs.")
	return pairs
}

// extendTrajectory tries to extend a given trajectory (currentT) with a diagnosis (d). It returns a map which maps all
// patients that follow the extended trajectory onto an index in their diagnosis lists.
func extendTrajectory(currentT *Trajectory, d int, minTime, maxTime float64) map[*Patient]int {
	result := map[*Patient]int{}
	for p, idx := range currentT.TrajMap {
		idx2 := countPatientTrajectory(p, idx, d, minTime, maxTime)
		if idx2 != -1 {
			result[p] = idx2
		}
	}
	return result
}

// BuildTrajectories calculates the trajectories for an experiment. The trajectories are constrained by: a
// minimum number of patients in the trajectory (minPatients), a maximum number of diagnoses in the trajectory (maxLength),
// a minimum number of diagnoses in the trajectory (minLength), a minimum RR for each diagnosis transition (minRR), and
// a list of filters.
func (exp *Experiment) BuildTrajectories(minPatients, maxLength, minLength int, minTime, maxTime, minRR float64,
	filters []TrajectoryFilter) []*Trajectory {
	fmt.Println("Building patient trajectories...")
	pairs := exp.selectDiagnosisPairs(minPatients, minRR)
	exp.Pairs = pairs
	var trajectories []*Trajectory
	var stack []*Trajectory
	for _, pair := range pairs {
		t := &Trajectory{
			Diagnoses:      []int{pair.First, pair.Second},
			PatientNumbers: []int{len(exp.DxDPatients[pair.First][pair.Second])},
			Patients:       [][]*Patient{exp.DxDPatients[pair.First][pair.Second]},
			TrajMap:        map[*Patient]int{},
		}
		for _, p := range exp.DxDPatients[pair.First][pair.Second] {
			_, idx := countPatientDiagnosisPair(p, pair.First, pair.Second, minTime, maxTime)
			t.TrajMap[p] = idx
		}
		stack = append(stack, t)
	}
	// divide the work
	result := parallel.RangeReduce(0, len(stack), 0, func(low, high int) interface{} {
		lstack := stack[low:high]
		var ltrajectories []*Trajectory
		tCtr := 0
		for {
			if len(lstack) == 0 {
				break
			}
			currentT := lstack[0]
			lstack = lstack[1:]
			// find potential extensions
			lastT := currentT.Diagnoses[len(currentT.Diagnoses)-1]
			ctr := 0
			for _, pair := range pairs {
				if pair.First == lastT && len(exp.DxDPatients[lastT][pair.Second]) >= minPatients {
					//patients := intersectPatients(currentT.Patients[len(currentT.Patients)-1], exp.DxDPatients[lastT][pair.Second])
					extendedTrajMap := extendTrajectory(currentT, pair.Second, minTime, maxTime)
					if len(extendedTrajMap) > minPatients {
						currentT.TrajMap = extendedTrajMap
						diagnoses := make([]int, len(currentT.Diagnoses))
						copy(diagnoses, currentT.Diagnoses)
						patientNumbers := make([]int, len(currentT.PatientNumbers))
						copy(patientNumbers, currentT.PatientNumbers)
						ps := make([][]*Patient, len(currentT.Patients))
						copy(ps, currentT.Patients)
						var patients []*Patient
						for p := range extendedTrajMap {
							patients = append(patients, p)
						}
						newT := &Trajectory{
							Diagnoses:      append(diagnoses, pair.Second), // should copy slice, could be updated many times...
							PatientNumbers: append(patientNumbers, len(patients)),
							Patients:       append(ps, patients),
						}
						// check if trajectory is finalized
						if len(newT.Diagnoses) >= maxLength {
							//newT.Patients = nil // help gc
							ltrajectories = append(ltrajectories, newT)
							tCtr++
						} else {
							ctr++
							lstack = append(lstack, newT)
						}
					}
				}
			}
			if ctr == 0 && len(currentT.Diagnoses) >= minLength { // no extension, finalize this trajectory
				ltrajectories = append(ltrajectories, currentT)
				tCtr++
			}
		}
		return ltrajectories
	}, func(result1, result2 interface{}) interface{} {
		r1 := result1.([]*Trajectory)
		r2 := result2.([]*Trajectory)
		for _, t := range r2 {
			r1 = append(r1, t)
		}
		return r1
	})
	trajectories = result.([]*Trajectory)
	fmt.Println("Found ", len(trajectories), " trajectories.")
	var filteredTrajectories []*Trajectory
	for idx, traj := range trajectories {
		keep := true
		for _, filter := range filters {
			if !filter(traj) {
				keep = false
				break
			}
		}
		if keep {
			traj.ID = idx
			filteredTrajectories = append(filteredTrajectories, traj)
		}
	}
	fmt.Println("Filtered down from: ", len(trajectories), " trajectories down to: ", len(filteredTrajectories),
		" trajectories.")
	exp.Trajectories = filteredTrajectories
	return filteredTrajectories
}
