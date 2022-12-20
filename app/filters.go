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
	"ptra/trajectory"
	"strings"
)

// cancerStageAggregator filters a set of patients to only include those that satisfy a given predicate that is applied
// on the patient's tumor information (which encodes cancer stages etc).
func cancerStageAggregator(predicate func(tInfo *TumorInfo) bool, tInfoMap map[string][]*TumorInfo) trajectory.PatientFilter {
	return func(p *trajectory.Patient) bool {
		//multiple tumor info entries per patient possible
		if tInfos, ok := tInfoMap[p.PIDString]; ok {
			if !ok {
				return false
			}
			tInfoToUseIndex := -1
			for i, tInfo := range tInfos { //go over all infos to grab the latest cancer stage that satisfies the predicate
				if predicate(tInfo) {
					tInfoToUseIndex = i
				}
			}
			if tInfoToUseIndex != -1 {
				// have a patient with specific cancer stage diagnosis
				// filter out diagnoses at later dates if possibly followed by other cancer stage
				if tInfoToUseIndex+1 < len(tInfos) {
					nextStageDate := tInfos[tInfoToUseIndex+1].Date
					newD := []*trajectory.Diagnosis{}
					for _, d := range p.Diagnoses {
						if trajectory.DiagnosisDateSmallerThan(d.Date, nextStageDate) {
							newD = append(newD, d)
						} else {
							continue
						}
					}
					p.Diagnoses = newD
				}
				return true
			}
		}
		return false
	}
}

// NMIBCAggregator checks all patients if they match the cancer criteria to be defined as non muscle invasive bladder
// cancer patients.
func NMIBCAggregator(tinfoMap map[string][]*TumorInfo) trajectory.PatientFilter {
	return cancerStageAggregator(func(tInfo *TumorInfo) bool {
		if tInfo.TStage == "Tis" || tInfo.TStage == "Ta" ||
			(tInfo.TStage == "T1" && tInfo.NStage == "N0" && tInfo.MStage == "M0") {
			return true
		}
		return false
	}, tinfoMap)
}

// MIBCAggregator checks all patients if they match the cancer criteria to be defined as muscle invasive bladder cancer
// patients.
func MIBCAggregator(tinfoMap map[string][]*TumorInfo) trajectory.PatientFilter {
	return cancerStageAggregator(func(tInfo *TumorInfo) bool {
		if tInfo.TStage == "T2" || tInfo.TStage == "T3" ||
			(tInfo.TStage == "T4" && tInfo.MStage == "M0" &&
				(tInfo.NStage == "N0" || tInfo.NStage == "N1" || tInfo.NStage == "N2" || tInfo.NStage == "N3")) {
			return true
		}
		return false
	}, tinfoMap)
}

// MUCAggregator checks all patients if they match the cancer criteria to be defined as metastisized bladder cancer
// patients.
func MUCAggregator(tinfoMap map[string][]*TumorInfo) trajectory.PatientFilter {
	return cancerStageAggregator(func(tInfo *TumorInfo) bool {
		if tInfo.MStage == "M0" {
			return true
		}
		return false
	}, tinfoMap)
}

// TaStageAggregator collects patients with stage Ta bladder cancer.
func TaStageAggregator(tinfoMap map[string][]*TumorInfo) trajectory.PatientFilter {
	return cancerStageAggregator(func(tInfo *TumorInfo) bool {
		if tInfo.TStage == "Ta" {
			return true
		}
		return false
	}, tinfoMap)
}

// T1StageAggregator collects patients with stage T1 bladder cancer.
func T1StageAggregator(tinfoMap map[string][]*TumorInfo) trajectory.PatientFilter {
	return cancerStageAggregator(func(tInfo *TumorInfo) bool {
		if tInfo.TStage == "T1" || tInfo.TStage == "T1a" || tInfo.TStage == "T1c" {
			return true
		}
		return false
	}, tinfoMap)
}

// TisStageAggregator collects patients with stage Tis bladder cancer.
func TisStageAggregator(tinfoMap map[string][]*TumorInfo) trajectory.PatientFilter {
	return cancerStageAggregator(func(tInfo *TumorInfo) bool {
		if tInfo.TStage == "Tis" {
			return true
		}
		return false
	}, tinfoMap)
}

// T2StageAggregator collects patients with stage T2 bladder cancer.
func T2StageAggregator(tinfoMap map[string][]*TumorInfo) trajectory.PatientFilter {
	return cancerStageAggregator(func(tInfo *TumorInfo) bool {
		if tInfo.TStage == "T2" || tInfo.TStage == "T2a" || tInfo.TStage == "T2b" || tInfo.TStage == "T2c" {
			return true
		}
		return false
	}, tinfoMap)
}

// T3StageAggregator collects patients with stage T3 bladder cancer.
func T3StageAggregator(tinfoMap map[string][]*TumorInfo) trajectory.PatientFilter {
	return cancerStageAggregator(func(tInfo *TumorInfo) bool {
		if tInfo.TStage == "T3" || tInfo.TStage == "T3a" || tInfo.TStage == "T3b" {
			return true
		}
		return false
	}, tinfoMap)
}

// T4StageAggregator collects patients with stage T4 bladder cancer.
func T4StageAggregator(tinfoMap map[string][]*TumorInfo) trajectory.PatientFilter {
	return cancerStageAggregator(func(tInfo *TumorInfo) bool {
		if tInfo.TStage == "T4" || tInfo.TStage == "T4a" || tInfo.TStage == "T4b" {
			return true
		}
		return false
	}, tinfoMap)
}

// N0StageAggregator collects patients with stage N0 bladder cancer.
func N0StageAggregator(tinfoMap map[string][]*TumorInfo) trajectory.PatientFilter {
	return cancerStageAggregator(func(tInfo *TumorInfo) bool {
		if tInfo.NStage == "N0" {
			return true
		}
		return false
	}, tinfoMap)
}

// N1StageAggregator collects patients with stage N1 bladder cancer.
func N1StageAggregator(tinfoMap map[string][]*TumorInfo) trajectory.PatientFilter {
	return cancerStageAggregator(func(tInfo *TumorInfo) bool {
		if tInfo.NStage == "N1" {
			return true
		}
		return false
	}, tinfoMap)
}

// N2StageAggregator collects patients with stage N2 bladder cancer.
func N2StageAggregator(tinfoMap map[string][]*TumorInfo) trajectory.PatientFilter {
	return cancerStageAggregator(func(tInfo *TumorInfo) bool {
		if tInfo.NStage == "N2" {
			return true
		}
		return false
	}, tinfoMap)
}

// N3StageAggregator collects patients with stage N0 bladder cancer.
func N3StageAggregator(tinfoMap map[string][]*TumorInfo) trajectory.PatientFilter {
	return cancerStageAggregator(func(tInfo *TumorInfo) bool {
		if tInfo.NStage == "N3" {
			return true
		}
		return false
	}, tinfoMap)
}

// M0StageAggregator collects patients with stage M0 bladder cancer.
func M0StageAggregator(tinfoMap map[string][]*TumorInfo) trajectory.PatientFilter {
	return cancerStageAggregator(func(tInfo *TumorInfo) bool {
		if tInfo.MStage == "M0" {
			return true
		}
		return false
	}, tinfoMap)
}

// M1StageAggregator collects patients with stage M1 bladder cancer.
func M1StageAggregator(tinfoMap map[string][]*TumorInfo) trajectory.PatientFilter {
	return cancerStageAggregator(func(tInfo *TumorInfo) bool {
		if tInfo.MStage == "M1" || tInfo.MStage == "M1a" || tInfo.MStage == "M1b" {
			return true
		}
		return false
	}, tinfoMap)
}

// CancerTrajectoryFilter filters trajectories down to trajectories that contain at least one diagnosis code that is
// considered to be cancer-related, e.g. containing the word "neoplasm".
func CancerTrajectoryFilter(exp *trajectory.Experiment) trajectory.TrajectoryFilter {
	//Determine all diagnosis codes that are cancer-related
	CancerRelatedMap := map[int]bool{}
	for did, medName := range exp.NameMap {
		medWords := strings.Split(medName, " ")
		cancerRelated := false
		for _, word := range medWords {
			if word == "neoplasm" || word == "Neoplasm" {
				cancerRelated = true
				break
			}
		}
		CancerRelatedMap[did] = cancerRelated
	}
	return func(t *trajectory.Trajectory) bool {
		for _, did := range t.Diagnoses {
			if CancerRelatedMap[did] {
				return true
			}
		}
		return false
	}
}

// BladderCancerTrajectoryFilter filters trajectories down to trajectorories with at least one diagnosis that is related
// to bladder cancer specifically, cf. ICD10 categories C67,C77,C78,C79 or procedures such as MVAC chemo, IVT treatment,
// or radical cystectomy.
func BladderCancerTrajectoryFilter(exp *trajectory.Experiment) trajectory.TrajectoryFilter {
	//Determine all internal diagnosis codes that are bladder cancer-related
	bladderCancerRelatedMap := map[int]bool{}
	for did, _ := range exp.NameMap {
		icdCode := exp.IdMap[did]
		if len(icdCode) >= 3 {
			subCode := icdCode[0:3]
			if subCode == "C67" || subCode == "C77" || subCode == "C78" || subCode == "C79" || // also check self-defined codes for treatments
				subCode == "C98" || subCode == "C99" || (len(icdCode) >= 4 && icdCode[0:4] == "C100") {
				bladderCancerRelatedMap[did] = true
			} else {
				bladderCancerRelatedMap[did] = false
			}
		}
	}
	return func(t *trajectory.Trajectory) bool {
		for _, did := range t.Diagnoses {
			if bladderCancerRelatedMap[did] {
				return true
			}
		}
		return false
	}
}
