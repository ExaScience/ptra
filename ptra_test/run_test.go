package ptra_test

import (
	"log/slog"
	"ptra/app"
	"ptra/cluster"
	"ptra/trajectory"
	"runtime"
	"testing"
)

//This test is the equivalent of running the command:
//./ptra fake-patients-simplified.csv icd10cm_tabular_2022.xml fake-diagnoses-simplified.csv
//		 ./exp1-test2/ --nofAgeGroups 6 --lvl 2 --maxYears 40.0 --minYears 0.5
//		 --minPatients 10 --maxTrajectoryLength 5 --minTrajectoryLength 3
//		 --name exp1-test2 --iter 100 --cluster --mclPath ~/tools/mcl/

func testMainLoop(t *testing.T) {
	slog.Info("***********************")
	slog.Info("** Test 2: Main loop **")
	slog.Info("***********************")
	//inputs
	pfile := "./data/fake-patients-simplified.csv"  //patient info
	hfile := "./data/icd10cm_tabular_2022.xml"      //icd 10 hierarchy
	dfile := "./data/fake-diagnoses-simplified.csv" //diagnosis files
	mclPath := "~/tools/mcl"
	afile := ""
	pfilters := []trajectory.PatientFilter{}
	eoid := []string{}
	//parameters
	n := 6   //age groups
	lvl := 2 //icd level
	maxYears := 40.0
	minYears := 0.5
	minPatients := 10
	maxTrajectoryLength := 5
	minTrajectoryLength := 3
	name := "exp1-test2"
	iter := 100
	//output
	output := "./exp1-test2/"
	slog.Info("Parsing inputs")
	//1. parse inputs into experiment object
	exp, _ := app.ParseData(name, pfile, dfile, hfile, afile, n, lvl, pfilters, eoid)
	//2. calculate relative risk ratios
	slog.Info("Calculating relative risk ratios")
	trajectory.InitializeExperimentRelativeRiskRatios(exp, minYears, maxYears, iter)
	//3. Build the trajectories
	slog.Info("Build the trajectories")
	rr := 1.0 //the minimum RR score to consider a diagnosis pair for trajectory building
	tfilters := []trajectory.TrajectoryFilter{}
	trajectory.BuildTrajectories(exp, minPatients, maxTrajectoryLength, minTrajectoryLength, minYears, maxYears, rr, tfilters)
	//4. Plot trajectories to file
	slog.Info("Plot trajectories to file")
	trajectory.PrintTrajectoriesToFile(exp, output)
	//5. Cluster when running on Linux
	if runtime.GOOS == "linux" {
		clusterGranularityList := []int{40, 60, 80, 100}
		slog.Info("MCL Clustering")
		cluster.ClusterTrajectoriesDirectly(exp, clusterGranularityList, output, mclPath)
	} else {
		slog.Info("Skipping clustering. Mcl not available on: " + runtime.GOOS)
	}
}
