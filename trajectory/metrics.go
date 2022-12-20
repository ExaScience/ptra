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

import "math"

// Collecting metrics for clusters of trajectories

// AgeAtDiagnosis calculates the age of a patient at a specific diagnosis
func AgeAtDiagnosis(p *Patient, DID int) int {
	yob := p.YOB
	var diagnosis *Diagnosis
	for _, d := range p.Diagnoses {
		if d.DID == DID {
			diagnosis = d
			break
		}
	}
	date := diagnosis.Date
	return date.Year - yob
}

// AgeAtEOI calculates the age of a patient at the event of interest (e.g. cancer diagnosis)
func AgeAtEOI(p *Patient) int {
	yob := p.YOB
	if p.EOIDate != nil {
		return p.EOIDate.Year - yob
	}
	return -1
}

// MetricsFromTrajectories computes:
// * mean age + standard deviation + median age for patients in the trajectories. Patients can occur in different
// trajectories. For mean age + sd + median, they will be counted as separate instances.
// * #patients per age category (normalized in percentages, not absolute numbers). Patients that occr in different
// trajectories will be counted as separate instances for these age categories.
// * #males, #females
// * mean survival time after event of interest
func MetricsFromTrajectories(trajectories []*Trajectory) (float64, float64, float64, float64, int, int) {
	meanAge := 0
	ctr := 0
	mCtr := 0
	fCtr := 0
	meanAgeOfEOI := 0
	ctr2 := 0
	for _, t := range trajectories {
		for _, p := range t.Patients[len(t.Patients)-1] { // patients in last diagnosis of the trajectory
			ctr++
			meanAge = meanAge + AgeAtDiagnosis(p, t.Diagnoses[len(t.Diagnoses)-1])
			if p.Sex == Male {
				mCtr++
			} else {
				fCtr++
			}
			ageEOI := AgeAtEOI(p)
			if ageEOI != -1 {
				meanAgeOfEOI = meanAgeOfEOI + ageEOI
				ctr2++
			}
		}
	}
	meanAgeF := float64(meanAge) / float64(ctr)
	meanAgeOfEOIF := float64(meanAgeOfEOI) / float64(ctr2)
	stdDev := 0.0
	stdDevEOI := 0.0
	for _, t := range trajectories {
		for _, p := range t.Patients[len(t.Patients)-1] { // patients in last diagnosis of the trajectory
			age := float64(AgeAtDiagnosis(p, t.Diagnoses[len(t.Diagnoses)-1]))
			stdDev = stdDev + ((meanAgeF - age) * (meanAgeF - age))
			ageEOI := float64(AgeAtEOI(p))
			if ageEOI != -1 {
				stdDevEOI = stdDevEOI + ((meanAgeOfEOIF - ageEOI) * (meanAgeOfEOIF - ageEOI))
			}
		}
	}
	stdDev = math.Sqrt(stdDev / float64(ctr))
	stdDevEOI = math.Sqrt(stdDevEOI / float64(ctr2))
	return meanAgeF, stdDev, meanAgeOfEOIF, stdDevEOI, mCtr, fCtr
}
