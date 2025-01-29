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
	"errors"
	"fmt"
	"github.com/imec-int/ptra/lib/utils"
	"os"
	"path"
	"runtime"
	"runtime/debug"
	"strconv"
	"strings"
)

type ExperimentParams struct {
	// required parameters
	Name             string
	PatientInfo      string // path to the file with patient information (ID, gender," + birthyear, etc)
	DiagnosisInfo    string // path to the file with diagnosis information (ID,descriptor, hierarchy, etc)
	PatientDiagnoses string // path to the file with patient diagnoses.
	OutputPath       string // path where output files are written to.

	// optional parameters
	NofAgeGroups         int
	Lvl                  int
	MaxYears             float64
	MinYears             float64
	MinPatients          int
	MaxTrajectoryLength  int
	MinTrajectoryLength  int
	ICD9ToICD10File      string
	Cluster              bool
	ClusterGranularities string
	Iter                 int
	RR                   float64
	SaveRR               string
	LoadRR               string
	PFilters             string
	TFilters             string
	TumorInfo            string
	TreatmentInfo        string
	NrOfThreads          int
}

// Run runs a TriNetX experiment with the given parameters.
func Run(args *ExperimentParams) (err error) {
	defer func() {
		// converts any panics into errors to avoid crashing the app
		if r := recover(); r != nil {
			fmt.Println("Recovered from panic during experiment: ", r)
			err = errors.New(fmt.Sprintf("%v", r))
			fmt.Println(string(debug.Stack()))
		}
	}()

	outputDir := path.Join(args.OutputPath, args.Name)
	err = os.MkdirAll(outputDir, 0700)
	if err != nil {
		return err
	}

	if args.NrOfThreads > 0 {
		runtime.GOMAXPROCS(args.NrOfThreads)
	}

	// start execution
	// 1. Parse input into experiment
	tinfo := map[string][]*TumorInfo{}
	if args.TumorInfo != "" {
		tinfo = ParsetTriNetXTumorData(args.TumorInfo) // need parsed patients to be able to parse tumor data file
	}

	exp, patients := ParseTriNetXData(args.Name, args.PatientInfo, args.PatientDiagnoses, args.DiagnosisInfo,
		args.TreatmentInfo, args.NofAgeGroups, args.Lvl, args.MinYears, args.MaxYears, args.ICD9ToICD10File, GetPatientFilters(args.PFilters, tinfo))

	// 2. Initialise relative risk ratios or load them from file from a previous run
	if args.LoadRR != "" {
		exp.LoadRRMatrix(args.LoadRR)
		exp.LoadDxDPatients(patients, fmt.Sprintf("%s.patients.csv", args.LoadRR))
	} else {
		exp.InitRR(args.MinYears, args.MaxYears, args.Iter)
	}
	if args.SaveRR != "" { // save RR matrix to file + DPatients
		exp.SaveRRMatrix(args.SaveRR)
		exp.SaveDxDPatients(fmt.Sprintf("%s.patients.csv", args.SaveRR))
	}

	// assist the gc and nil some exp data that is no longer needed after initializing RR
	exp.Cohorts = nil
	exp.DPatients = nil

	// 3. Build the trajectories
	exp.BuildTrajectories(args.MinPatients, args.MaxTrajectoryLength, args.MinTrajectoryLength, args.MinYears, args.MaxYears, args.RR,
		GetTrajectoryFilters(args.TFilters, exp))

	// 4. Plot trajectories to file
	exp.PrintTrajectoriesToFile(outputDir)
	fmt.Println("Collected trajectories: ")
	for i := 0; i < utils.MinInt(len(exp.Trajectories), 100); i++ {
		LogTrajectory(exp.Trajectories[i], exp)
	}

	// 5. Perform clustering
	if args.Cluster {
		var clusterGranularityList []int
		for _, g := range strings.Split(args.ClusterGranularities, ",") {
			gi, _ := strconv.ParseInt(g, 10, 0)
			clusterGranularityList = append(clusterGranularityList, int(gi))
		}
		clusteringErr := ClusterTrajectories(exp, clusterGranularityList, outputDir)
		if clusteringErr != nil {
			return clusteringErr
		}
	}

	return nil
}
