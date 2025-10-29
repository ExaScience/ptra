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

package main

import (
	"bytes"
	"log"
	"ptra/app"
	"ptra/cluster"
	"ptra/trajectory"
	"ptra/utils"
	"strconv"
	"strings"

	//"bytes"
	"flag"
	"fmt"
	"io/ioutil"
	"path/filepath"

	//"log"
	"os"
	"runtime"
)

/*
Ptra is a tool for patient trajectory analysis.

Usage:
	ptra pfile ifile dfile path [flags]

Example:
	ptra ICD10 patient.csv icd10cm_tabular_2022.xml diagnosis.csv ./MIBC_tfiltered/ --nofAgeGroups 10 --lvl 2
	--maxYears 5 --minYears 0.001 --minPatients 50 --maxTrajectoryLength 5 --minTrajectoryLength 3 --name MICB_tfiltered
	--ICD9ToICD10File ICD_9_to_10.json --iter 400 --RR 1.0 --saveRR MIBC_tfiltered.csv --cluster
	--mclPath /home/caherzee/tools/mcl/ --clusterGranularities 40,60,80,100

The flags are:

--nofAgeGroups nr
	The number of age groups to consider when dividing the population into cohorts based on age. The tool automatically
	detects the oldest and youngest patients from input. The number of age groups is used to divide the minumum birth
	year and maximum birth year into age ranges. E.g. min birth year: 1900, max birth year: 2020, and nr of age groups:
	10, will create 10 cohorts where age ranges from [0,12],[12,22],...[108,120].
--lvl nr
	Sets the ICD10 level for the diagnosis input [0-6]. The tool maps all ICD10 codes to a medical meaningful term based on
	this chosen level. ICD10 codes of lower levels may be combined into the same code of a higher level. E.g. A00.0
	Cholera due to Vibrio cholerae 01, biovar cholerae and A00.1 Cholera due to Vibrio cholerae 01, biovar eltor are lvl
	3 codes and may be collapsed to A00 Cholera in lvl 2, or A00-A09 Intestinal infectious diseases in lvl 1, or A00-B99
	Certain infectious and parasitic diseases in lvl 0.
--minPatients nr
	Sets the minimum required number of patients in a trajectory.
--maxYears nr
	Sets the maximum number of years between subsequent diagnoses to be considered for inclusion in a trajectory. E.g.
	0.5 for half a year.
--minYears nr
	Sets the minimum number of years between subsequent diagnoses to be considered for inclusion in a trajectory. E.g.
	0.5 for half a year.
--maxTrajectoryLength nr
	Sets the maximum length of trajectories to be included in the output. E.g. 5 for trajectories with maximum 5
	diagnoses.
--minTrajectoryLength nr
	Sets the minimum length of trajectories to be included in the output. E.g. 3 for trajectories with minimum 3
	diagnoses.
--name string
	Sets the name of the experiment. This name is used to generate names for output files.
--cluster
	If this flag is passed, the computed trajectories are clustered and the clusters are outputted to file.
--clusterGranularities nr, nr, .., nr
	The cluster granularities to try. This impacts how many clusters the algorithm tries to find. Default: 40,60,80,100
--mclPath
	Sets the path where the mcl binaries can be found.
--iter nr
	Sets the number of iterations to be used in the sampling experiments for calculating relative risk ratios. If iter
	is 400, the calculated p-values are within 0.05 of the true p-values. For iter = 10000, the true p-values are within
	0.01 of the true p-values. The higher the number of iterations, the higher the runtime.
--saveRR file
	Save the RR matrix, a matrix that represents the RR calculated from the population for each possible combination of
	ICD10 diagnosis pairs. This matrix can be loaded in other ptra runs to avoid recalculating the RR scores. This can
	be useful if parameters want to be explored that do not impact the RR calculation itself. Only iter, maxYears and
	minYears, and filters influence RR calculation. Variations of other parameters for constructing trajectories from RR
	scores, such as maxTrajectoryLenght, minTrajectoryLength, minPatients, RR etc might be explored in other runs.
--loadRR file
	Load the RR matrix from file. Such a file must be created by a previous run of ptra with the --saveRR flag.
--pfilters age+:nr | age-:nr | sex:male | sex:female
	A list of filters for selecting patients from which to derive trajectories. The form of a filter is tag:value. E.g.
	sex:male says to select only male patients.
--eoi icd10,icd10,...icd10
	A list of icd10 codes for identifying the event of interest. The event of interest can be used for analyzing or
	filtering trajectories
--tfilters cat:[neoplasm | bc], code:icd10;icd10;...icd10
	A list of filters for reducing the output of trajectories. The form to specify a filter is tag:value.
    E.g. cat:neoplasm only outputs trajectories where there is at least one diagnosis related to cancer. cat:bc only
	outputs trajectories where at least one diagnosis is related to bladder cancer. The code: tag can be used to list
	a number of ICD10 codes of which at least one should occur in the trajectory.
*/

const (
	programVersion = 0.11
	programName    = "ptra"
)

func programMessage() string {
	return fmt.Sprint(programName, " version ", programVersion, " compiled with ", runtime.Version())
}

const ptraHelp = "\nptra parameters:\n" +
	"ptra patientInfoFile diagnosisInfoFile diagnosesFile outputPath \n" +
	"[--nofAgeGroups nr]\n" +
	"[--lvl nr]\n" +
	"[--minPatients nr]\n" +
	"[--maxYears nr]\n" +
	"[--minYears nr]\n" +
	"[--maxTrajectoryLength nr]\n" +
	"[--minTrajectoryLength nr]\n" +
	"[--name string]\n" +
	"[--eoid icd10,icd10,...,icd10\n" +
	"[--cluster]\n" +
	"[--clusterGranularities nr,nr,...,nr]\n" +
	"[--mclPath string]\n" +
	"[--iter nr]\n" +
	"[--saveRR file]\n" +
	"[--loadRR file]\n" +
	"[--pfilters age+:nr,age-:nr,[sex:male | sex:female] ]\n" +
	"[--tfilters cat:[neoplasms | bc],code:icd10;...;icd10]\n" +
	"[--nrOfThreads nr]\n" +
	"[--actFile file]\n"

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

func getPatientFilter(s string) trajectory.PatientFilter {
	id := func(p *trajectory.Patient) bool { return true }
	switch s {
	case "id":
		return id
	case "EOI-":
		return trajectory.EOIAfterFilter()
	case "EOI+":
		return trajectory.EOIBeforeFilter()
	}
	fs := strings.Split(s, ":")
	if len(fs) > 1 {
		switch fs[0] {
		case "age+":
			{
				age, err := strconv.Atoi(fs[1])
				if err != nil {
					panic("Age passed to age+ not a well-formed integer: " + fs[1])
				}
				return trajectory.AgeAboveAggregator(age)
			}
		case "age-":
			{
				age, err := strconv.Atoi(fs[1])
				if err != nil {
					panic("Age passe to age- not a well-formed integer: " + fs[1])
				}
				return trajectory.AgeLessThanAggregator(age)
			}
		case "sex":
			{
				gender := fs[1]
				if gender == "male" {
					return trajectory.FemaleFilter()
				}
				if gender == "female" {
					return trajectory.MaleFilter()
				}
				panic("Unknown gender: " + gender)
			}
		default:
			panic("Unknown patient filter: " + s)
		}
	} else {
		switch s {
		case "eoi-":
			return trajectory.EOIAfterFilter()
		case "eoi+":
			return trajectory.EOIBeforeFilter()
		case "id":
			return id
		default:
			panic("Unknown patient filter: " + s)
		}
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

func getTrajectoryFilter(s string, exp *trajectory.Experiment) trajectory.TrajectoryFilter {
	id := func(t *trajectory.Trajectory) bool { return true }
	ss := strings.Split(s, ":")
	switch ss[0] {
	case "cat":
		switch ss[1] {
		case "neoplasm":
			return app.CancerTrajectoryFilter(exp)
		case "bc":
			return app.BladderCancerTrajectoryFilter(exp)
		default:
			panic("Unknown category for trajectory filter: cat:" + ss[1])
		}
	case "code":
		codes := strings.Split(ss[1], ";")
		return app.Icd10TrajectoryFilter(exp, codes)
	case "id":
		return id
	default:
		panic("Unknown trajectory filter: " + ss[0])
	}
}

func getTrajectoryFilters(f string, exp *trajectory.Experiment) []trajectory.TrajectoryFilter {
	fs := strings.Split(f, ",")
	result := []trajectory.TrajectoryFilter{}
	for _, f := range fs {
		result = append(result, getTrajectoryFilter(f, exp))
	}
	return result
}

func main() {
	var (
		// required parameters
		patientInfo      string //The file with patient information (ID, gender," + birthyear, etc)
		diagnosisInfo    string //The file with diagnosis information (ID,descriptor, hierarchy, etc)
		patientDiagnoses string //The file with patient diagnoses.
		outputPath       string //The path where output files are written.
		// optional flags
		nofAgeGroups         int
		lvl                  int
		maxYears             float64
		minYears             float64
		minPatients          int
		maxTrajectoryLength  int
		minTrajectoryLength  int
		name                 string
		clust                bool
		mclPath              string
		clusterGranularities string
		iter                 int
		rr                   float64
		saveRR               string
		loadRR               string
		pfilters             string
		tfilters             string
		tumorInfo            string
		treatmentInfo        string
		nrOfThreads          int
		eoid                 string
		actFile              string
	)
	var flags flag.FlagSet
	// options for the ptra command
	flags.IntVar(&nofAgeGroups, "nofAgeGroups", 6, "The population data is divided in cohorts in"+
		"terms of age groups to calculate relative risk ratios of diagnosis pairs. This parameters configures how"+
		"many age groups to use")
	flags.IntVar(&nrOfThreads, "nrOfThreads", 0, "The number of threads ptra uses.")
	flags.IntVar(&lvl, "lvl", 3, "Diagnosis codes are organised in a hierarchy of diagnosis "+
		"descriptors. The level says which descriptor in the hiearchy to use for trajectory building.")
	flags.Float64Var(&maxYears, "maxYears", 5.0, "The maximum number of years between diagnosis "+
		"A and B to consider the diagnosis pair A->B in a trajectory.")
	flags.Float64Var(&minYears, "minYears", 0.5, "The minimum number of years between diagnisis "+
		"A and B to consider the diagnosis pair A->B in a trajectory.")
	flags.IntVar(&minPatients, "minPatients", 1000, "The minimum number of patients for the last "+
		"diagnosis in a trajectory")
	flags.IntVar(&maxTrajectoryLength, "maxTrajectoryLength", 5, "The maximum number of diagnoses"+
		" in a trajectory")
	flags.IntVar(&minTrajectoryLength, "minTrajectoryLength", 3, "The minimum number of "+
		"diagnoses in a trajectory")
	flags.StringVar(&name, "name", "exp1", "The name of the run. This is used to generate the "+
		"names of the output files.")
	flags.StringVar(&eoid, "eoid", "", "The icd10 codes of the events of interest to track")
	flags.BoolVar(&clust, "cluster", false, "Cluster the trajectories using MCL and output "+
		"the results")
	flags.StringVar(&mclPath, "mclPath", "", "The path to the mcl binary.")
	flags.StringVar(&clusterGranularities, "clusterGranularities", "40,60,80,100", "The "+
		"granularities used for the mcl clustering step.") // recommended 14,20,40,60
	flags.IntVar(&iter, "iter", 10000, "The minimum number of sampling iterations "+
		"diagnosis in a trajectory")
	flags.Float64Var(&rr, "RR", 1.0, "The minimum RR score for considering pairs.")
	flags.StringVar(&saveRR, "saveRR", "", "Save the RR matrix to a file so it can be loaded for "+
		"later runs")
	flags.StringVar(&loadRR, "loadRR", "", "Load the RR matrix from a given file instead of "+
		"calculating it from scratch.")
	flags.StringVar(&pfilters, "pfilters", "id", "A list of pfilters to restrict analysis on specific "+
		"patients.")
	flags.StringVar(&tumorInfo, "tumorInfo", "", "A file with information about the tumor stages.")
	flags.StringVar(&treatmentInfo, "treatmentInfo", "", "A file with information about patient cancer stages.")
	flags.StringVar(&tfilters, "tfilters", "id", "A list of filters to restrict output of trajectories."+
		"A list of options of the form \"code:icd;icd;icd\" or \"cat:[neoplasm | bc]\"")
	flags.StringVar(&actFile, "actFile", "", "A files with ACT codes")
	// parse optional arguments
	parseFlags(flags, 5, ptraHelp)
	// parse required arguments
	patientInfo = getFileName(os.Args[1], ptraHelp)
	diagnosisInfo = getFileName(os.Args[2], ptraHelp)
	patientDiagnoses = getFileName(os.Args[3], ptraHelp)
	outputPath, _ = filepath.Abs(getFileName(os.Args[4], ptraHelp))
	outputPath = outputPath + string(filepath.Separator)
	fmt.Println("Output path: ", outputPath)
	// create output directory
	err := os.MkdirAll(filepath.Dir(outputPath), 0700)
	if err != nil {
		panic(err)
	}
	// build an output command line
	var command bytes.Buffer
	fmt.Fprint(&command, os.Args[0], " ", patientInfo, " ", diagnosisInfo, " ", patientDiagnoses,
		" ", outputPath)
	fmt.Fprint(&command, " --nofAgeGroups ", nofAgeGroups)
	fmt.Fprint(&command, " --lvl ", lvl)
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
	fmt.Fprint(&command, " --tfilters ", tfilters)
	if nrOfThreads > 0 {
		runtime.GOMAXPROCS(nrOfThreads)
		fmt.Fprint(&command, " --nrOfThreads ", nrOfThreads)
	}
	if eoid != "" {
		fmt.Fprint(&command, " --eoid ", eoid)
	}
	if actFile != "" {
		fmt.Fprint(&command, " --actFile ", actFile)
	}
	// start execution
	log.Println(programMessage())
	log.Println("Executing command:\n", command.String())
	//1. Parse inputs into experiment
	exp, patients := app.ParseData("exp1", patientInfo, patientDiagnoses, diagnosisInfo, actFile,
		nofAgeGroups, lvl, getPatientFilters(pfilters), strings.Split(eoid, ","))
	//2. Initialise relative risk ratios or load them from file from a previous run
	if loadRR != "" {
		trajectory.LoadRRMatrix(exp, loadRR)
		trajectory.LoadDxDPatients(exp, patients, fmt.Sprintf("%s.patients.csv", loadRR))
	} else {
		trajectory.InitializeExperimentRelativeRiskRatios(exp, minYears, maxYears, iter)
	}
	if saveRR != "" { //save RR matrix to file + DPatients
		trajectory.SaveRRMatrix(exp, saveRR)
		trajectory.SaveDxDPatients(exp, fmt.Sprintf("%s.patients.csv", saveRR))
	}
	// assist the gc and nil some exp data that is no longer needed after initializing RR
	exp.Cohorts = nil
	exp.DPatients = nil
	//3. Build the trajectories
	trajectory.BuildTrajectories(exp, minPatients, maxTrajectoryLength, minTrajectoryLength, minYears, maxYears, rr,
		getTrajectoryFilters(tfilters, exp))
	//4. Plot trajectories to file
	trajectory.PrintTrajectoriesToFile(exp, outputPath)
	fmt.Println("Collected trajectories: ")
	for i := 0; i < utils.MinInt(len(exp.Trajectories), 100); i++ {
		trajectory.PrintTrajectory(exp.Trajectories[i], exp)
	}
	//5. Perform clustering
	if clust {
		var clusterGranularityList []int
		for _, g := range strings.Split(clusterGranularities, ",") {
			gi, _ := strconv.ParseInt(g, 10, 0)
			clusterGranularityList = append(clusterGranularityList, int(gi))
		}
		fmt.Println("MCL Clustering:")
		//ClusterTrajectories(exp, clusterGranularityList, outputPath, mclPath)
		cluster.ClusterTrajectoriesDirectly(exp, clusterGranularityList, outputPath, mclPath)
	}
}
