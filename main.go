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
	--ICD9ToICD10File ICD_9_to_10.json --iter 400 --RR 1.0 --tumorInfo tumor.csv --saveRR MIBC_tfiltered.csv --cluster
	--mclPath /home/caherzee/tools/mcl/ --clusterGranularities 40,60,80,100 --pfilters "MIBC" --tumorInfo tumor.csv
	--tfilters "bc" --treatmentInfo treatments.csv

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
--ICD9ToICD10File file
	A json file that provides a mapping from ICD9 to ICD10 codes. The input may be mixed ICD9 and ICD10 codes. With this
	mapping, the tool can automatically convert all diagnosis codes to ICD10 codes for analysis.
--cluster
	If this flag is passed, the computed trajectories are clustered and the clusters are outputted to file.
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
--pfilters age70+ | age70- | male | female | Ta | T0 | Tis | T1 | T2 | T3 | T4 | N0 | N1 | N2 | N3 | M0 | M1 |NMIBC | MIBC | mUC
	A list of filters for selecting patients from whitch to derive trajectories.
--tumorInfo file
	A file with information about patients and their tumors. This file contains annotations about the stage of the
	bladder cancer at a specific time. Cf. TriNetX tumor table. This information is used by filters.
--tfilters neoplasm | bc
	A list of filters for reducing the output of trajectories. E.g. neoplasm only outputs trajectories where there is at
	least one diagnosis related to cancer. bc only outputs trajectories where one diagnosis is (assuming) related to
	bladder cancer.
--treatmentInfo file
	A file with information about patients and their treatments, e.g. MVAC,radical cystectomy, etc. If this file is
	passed, the treatments will be used as diagnostic codes to calculated trajectories.
*/

const (
	programVersion = 0.1
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
	"[--ICD9ToICD10File file]\n" +
	"[--cluster]\n" +
	"[--mclPath string]\n" +
	"[--iter nr]\n" +
	"[--saveRR file]\n" +
	"[--loadRR file]\n" +
	"[--pfilters age70+ | age70- | male | female | Ta | T0 | Tis | T1 | T2 | T3 | T4 | N0 | N1 | N2 | N3 | M0 | M1 |" +
	"NMIBC | MIBC | mUC ]\n" +
	"[--tumorInfo file]\n" +
	"[--tfilters neoplasm | bc]\n" +
	"[--treatmentInfo file]\n" +
	"[--nrOfThreads nr]\n"

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

func getPatientFilter(s string, tinfo map[string][]*app.TumorInfo) trajectory.PatientFilter {
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
	case "Ta":
		return app.TaStageAggregator(tinfo)
	case "T1":
		return app.T1StageAggregator(tinfo)
	case "Tis":
		return app.TisStageAggregator(tinfo)
	case "T2":
		return app.T2StageAggregator(tinfo)
	case "T3":
		return app.T3StageAggregator(tinfo)
	case "T4":
		return app.T4StageAggregator(tinfo)
	case "N0":
		return app.N0StageAggregator(tinfo)
	case "N1":
		return app.N1StageAggregator(tinfo)
	case "N2":
		return app.N2StageAggregator(tinfo)
	case "N3":
		return app.N3StageAggregator(tinfo)
	case "M0":
		return app.M0StageAggregator(tinfo)
	case "M1":
		return app.M1StageAggregator(tinfo)
	case "EOI-":
		return trajectory.EOIAfterFilter()
	case "EOI+":
		return trajectory.EOIBeforeFilter()
	case "MIBC":
		return app.MIBCAggregator(tinfo)
	case "NMIBC":
		return app.NMIBCAggregator(tinfo)
	case "mUC":
		return app.MUCAggregator(tinfo)
	default:
		return id
	}
}

func getPatientFilters(f string, tinfo map[string][]*app.TumorInfo) []trajectory.PatientFilter {
	fs := strings.Split(f, ",")
	result := []trajectory.PatientFilter{}
	for _, f := range fs {
		result = append(result, getPatientFilter(f, tinfo))
	}
	return result
}

func getTrajectoryFilter(s string, exp *trajectory.Experiment) trajectory.TrajectoryFilter {
	id := func(t *trajectory.Trajectory) bool { return true }
	switch s {
	case "neoplasm":
		return app.CancerTrajectoryFilter(exp)
	case "bc":
		return app.BladderCancerTrajectoryFilter(exp)
	default:
		return id
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
		ICD9ToICD10File      string
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
	flags.StringVar(&ICD9ToICD10File, "ICD9ToICD10File", "", "A json file that maps ICD9 to "+
		"ICD10 codes.")
	flags.BoolVar(&clust, "cluster", false, "Cluster the trajectories using MCL and output "+
		"the results")
	flags.StringVar(&mclPath, "mclPath", "/usr/bin/mcl", "The path to the mcl binary.")
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
	flags.StringVar(&tfilters, "tfilters", "id", "A list of pfilters to restrict output of trajectories")
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
	fmt.Fprint(&command, " --ICD9ToICD10File ", ICD9ToICD10File)
	fmt.Fprint(&command, " --iter ", iter)
	fmt.Fprint(&command, " --RR ", rr)
	fmt.Fprint(&command, " --tumorInfo ", tumorInfo)
	fmt.Fprint(&command, " --treatmentInfo ", treatmentInfo)
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
	// start execution
	log.Println(programMessage())
	log.Println("Executing command:\n", command.String())
	//1. Parse inputs into experiment
	// Parse Tumor info
	tinfo := map[string][]*app.TumorInfo{} // filterInfo is a variable to pass around filter-specific information. E.g. parsed tumor data for the tumor stage filter.
	if tumorInfo != "" {
		tinfo = app.ParsetTriNetXTumorData(tumorInfo) // need parsed patients to be able to parse tumor data file
	}
	exp, patients := app.ParseTriNetXData("exp1", patientInfo, patientDiagnoses, diagnosisInfo,
		treatmentInfo, nofAgeGroups, lvl, minYears, maxYears, ICD9ToICD10File, getPatientFilters(pfilters, tinfo))
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
