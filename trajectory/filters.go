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

package trajectory

// PatientFilter prescribes a function type for implementing filters on TriNetX patients, to be able to calculate
// trajectories for specific cohorts. E.g. male patients, patients <70 years, patients with specific cancer stage, etc.
type PatientFilter func(patient *Patient) bool

// TrajectoryFilter is a type to define a trajectory filter function. Such filters take as input a trajectory and must
// return a bool as output that determines if a trajectory passes a filter or not.
type TrajectoryFilter func(t *Trajectory) bool

func ApplyPatientFilter(filter PatientFilter, pMap *PatientMap) *PatientMap {
	newPMap := &PatientMap{PIDStringMap: map[string]int{}, PIDMap: map[int]*Patient{}, Ctr: pMap.Ctr}
	for pid, p := range pMap.PIDMap {
		if filter(p) {
			newPMap.PIDStringMap[p.PIDString] = pid
			newPMap.PIDMap[pid] = p
			if p.Sex == Male {
				newPMap.MaleCtr++
			} else {
				newPMap.FemaleCtr++
			}
		}
	}
	return newPMap
}

func ApplyPatientFilters(filters []PatientFilter, pMap *PatientMap) *PatientMap {
	newPMap := &PatientMap{PIDStringMap: map[string]int{}, PIDMap: map[int]*Patient{}, Ctr: pMap.Ctr}
	for pid, p := range pMap.PIDMap {
		res := true
		for _, filter := range filters {
			res = filter(p) && res
			if !res {
				break
			}
		}
		if res {
			newPMap.PIDStringMap[p.PIDString] = pid
			newPMap.PIDMap[pid] = p
			if p.Sex == Male {
				newPMap.MaleCtr++
			} else {
				newPMap.FemaleCtr++
			}
		}
	}
	return newPMap
}

// SexFilter removes all patients of the given sex.
func SexFilter(sex int) PatientFilter {
	return func(p *Patient) bool {
		return p.Sex != sex
	}
}

// MaleFilter removes all male patients.
func MaleFilter() PatientFilter {
	return SexFilter(Male)
}

// FemaleFilter removes all female patients.
func FemaleFilter() PatientFilter {
	return SexFilter(Female)
}

// EOIFilter removes all diagnoses for patients that satisfy a given predicate
func EOIFilter(test func(d1, d2 DiagnosisDate) bool) PatientFilter {
	return func(p *Patient) bool {
		if p.EOIDate == nil { //skip patients without EOIDate
			return false
		}
		newD := []*Diagnosis{}
		for _, d := range p.Diagnoses {
			if test(d.Date, *p.EOIDate) {
				break
			}
			newD = append(newD, d)
		}
		p.Diagnoses = newD
		if newD == nil {
			return false
		}
		return true
	}
}

// EOIBeforeFilter removes all diagnoses before the event of interest date
func EOIBeforeFilter() PatientFilter {
	return EOIFilter(func(d1, d2 DiagnosisDate) bool { return DiagnosisDateSmallerThan(d1, d2) })
}

// EOIAfterFilter removes all diagnoses after the event of interest date
func EOIAfterFilter() PatientFilter {
	return EOIFilter(func(d1, d2 DiagnosisDate) bool { return DiagnosisDateSmallerThan(d2, d1) })
}

// ageLessAggregator collects all patients younger than a specific age or trims down their data up until that age.
func ageLessAggregator(age int) PatientFilter {
	return func(p *Patient) bool {
		fYear := p.YOB + age - 1 // last year with diagnosis accepted
		//remove all diagnoses past a specific age
		newD := []*Diagnosis{}
		for _, d := range p.Diagnoses {
			if d.Date.Year > fYear {
				break
			}
			newD = append(newD, d)
		}
		p.Diagnoses = newD
		if len(newD) == 0 {
			return false
		}
		return true
	}
}

// ageAboveAggretator collects all patients older than a specific age and removes all diagnoses before that date.
func ageAboveAggregator(age int) PatientFilter {
	return func(p *Patient) bool {
		mYear := p.YOB + age // min year with diagnosis accepted
		//remove all diagnoses before a specific age
		newD := []*Diagnosis{}
		for _, d := range p.Diagnoses {
			if d.Date.Year <= mYear {
				continue
			}
			newD = append(newD, d)
		}
		p.Diagnoses = newD
		if len(newD) == 0 {
			return false
		}
		return true
	}
}

// LessThanSeventyAggregator collects all patients below a specific age.
func LessThanSeventyAggregator() PatientFilter {
	return ageLessAggregator(70)
}

// AboveSeventyAggregator collects all patients above a specific age.
func AboveSeventyAggregator() PatientFilter {
	return ageAboveAggregator(70)
}
