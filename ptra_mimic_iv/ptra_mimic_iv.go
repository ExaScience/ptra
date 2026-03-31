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
	"flag"
	"fmt"
	"io/ioutil"
	"log"
	"os"
	"path/filepath"
	"ptra/app"
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
	"[--ICD9ToICD10File file]          JSON or CSV mapping ICD-9 to ICD-10 (for diagnoses_icd)\n" +
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
		checkICD9            bool
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
	flags.StringVar(&icd9ToIcd10File, "ICD9ToICD10File", "", "JSON or CSV file mapping ICD-9 to ICD-10 codes (for diagnoses_icd source).")
	flags.BoolVar(&checkICD9, "checkICD9", false, "Check ICD-9 conversion coverage against diagnoses_icd.csv and exit.")

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

	// --checkICD9: run coverage check and exit
	if checkICD9 {
		dICDFile := filepath.Join(hospDir, "diagnoses_icd.csv")
		if _, err := os.Stat(dICDFile); os.IsNotExist(err) {
			fmt.Fprintf(os.Stderr, "Error: diagnoses_icd.csv not found in %s\n", hospDir)
			os.Exit(1)
		}
		if icd9ToIcd10File == "" {
			fmt.Fprintf(os.Stderr, "Error: --checkICD9 requires --ICD9ToICD10File\n")
			os.Exit(1)
		}
		app.CheckICD9Coverage(dICDFile, icd9ToIcd10File)
		return
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
		exp, patients = app.ParseHCPCSEventsWithPatients(hcpcsFile, patientsFile, nofAgeGroups,
			getPatientFilters(pfilters))
	case "diagnoses_icd":
		exp, patients = app.ParseDiagnosesICDWithPatients(diagnosesICDFile, patientsFile,
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
