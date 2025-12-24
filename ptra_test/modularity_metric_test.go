package ptra

import (
	"log/slog"
	"ptra/app"
	"ptra/cluster"
	"ptra/trajectory"
	"testing"
)

func createFakeDataExp() *trajectory.Experiment {
	//inputs
	pfile := "./data/fake-patients-simplified.csv"  //patient info
	hfile := "./data/icd10cm_tabular_2022.xml"      //icd 10 hierarchy
	dfile := "./data/fake-diagnoses-simplified.csv" //diagnosis files
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
	//1. parse inputs into experiment object
	exp, _ := app.ParseData(name, pfile, dfile, hfile, afile, n, lvl, pfilters, eoid)
	//2. calculate relative risk ratios
	trajectory.InitializeExperimentRelativeRiskRatios(exp, minYears, maxYears, iter)
	//3. Build the trajectories
	rr := 1.0 //the minimum RR score to consider a diagnosis pair for trajectory building
	tfilters := []trajectory.TrajectoryFilter{}
	trajectory.BuildTrajectories(exp, minPatients, maxTrajectoryLength, minTrajectoryLength, minYears, maxYears, rr, tfilters)
	return exp
}

func TestParseMCLFile(t *testing.T) {
	slog.Info("********************************************")
	slog.Info("** Test 3: Extended modularity clustering **")
	slog.Info("********************************************")
	file := "./data/fake-cluster-mcl.fake-data.txt"
	exp := createFakeDataExp()
	graph := cluster.ParseMCLGraph(exp, file)
	m := graph.ExtendedModularityMetric()
	slog.Info("Computed extended modularity metric: ")
	slog.Info("Value", "floatKey", m)
}
