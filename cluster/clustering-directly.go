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

package cluster

import (
	"bytes"
	"encoding/csv"
	"fmt"
	"io"
	"log"
	"log/slog"
	"os"
	"os/exec"
	"path/filepath"
	"ptra/trajectory"
	"ptra/utils"
	"strconv"
)

// jaccardTrajectory computes the Jaccard similarity coefficient for two given trajectories.
func jaccardTrajectory(t1, t2 *trajectory.Trajectory) float64 {
	// intersect t1 and t2
	n := 0
	for _, d1 := range t1.Diagnoses {
		if utils.MemberInt(d1, t2.Diagnoses) {
			n++
		}
	}
	nt1 := len(t1.Diagnoses)
	nt2 := len(t2.Diagnoses)
	return float64(n) / (float64(nt1) + float64(nt2) - float64(n))
}

// SzymkiewiczSimpsonTrajectory computes the Szymkiewicz-Simpson similarity coefficient for two given trajectories.
func SzymkiewiczSimpsonTrajectory(t1, t2 *trajectory.Trajectory) float64 {
	n := 0
	for _, d1 := range t1.Diagnoses {
		if utils.MemberInt(d1, t2.Diagnoses) {
			n++
		}
	}
	nt1 := len(t1.Diagnoses)
	nt2 := len(t2.Diagnoses)
	return float64(n) / float64(utils.MinInt(nt1, nt2))
}

// SorensenDiceTrajectory computes the SorensenDice similarity coefficient for two given trajectories.
func SorensenDiceTrajectory(t1, t2 *trajectory.Trajectory) float64 {
	n := 0
	for _, d1 := range t1.Diagnoses {
		if utils.MemberInt(d1, t2.Diagnoses) {
			n++
		}
	}
	nt1 := len(t1.Diagnoses)
	nt2 := len(t2.Diagnoses)
	return float64(2*n) / (float64(nt1 + nt2))
}

// convertTrajectoriesToAbcFormat compute the jaccard between each trajectory and writes out the result to file.
// Streaming algorithm to avoid pressure on memory.
func convertTrajectoriesToAbcFormat(exp *trajectory.Experiment, name string) {
	//create output file
	file, err := os.Create(name)
	if err != nil {
		log.Panic(err)
	}
	defer func() {
		if err := file.Close(); err != nil {
			log.Panic(err)
		}
	}()
	// compute the jacard index for the trajectories
	for i, t1 := range exp.Trajectories {
		t1.ID = i
		for j := i + 1; j < len(exp.Trajectories); j++ {
			t2 := exp.Trajectories[j]
			t2.ID = j
			coeff := jaccardTrajectory(t1, t2)
			fmt.Fprintf(file, "%d\t%d\t%f\n", i, j, coeff)
		}
	}
}

// ClusterTrajectoriesDirectly performs clustering of the trajectories that have been calculated for a given experiment.
// It does a pairwise comparison of all trajectories by calculating the jaccard similarity coefficients. Subsequently,
// MCL clustering is used to group the trajectories by jaccard similarity into clusters.
func ClusterTrajectoriesDirectly(exp *trajectory.Experiment, granularities []int, path, pathToMcl string) {
	slog.Info("Clustering trajectories directly with MCL")
	// convert trajectories to abc format for the mcl tool
	dirName := fmt.Sprintf("%s-clustered-trajectories/", exp.Name)
	workingDir := filepath.Join(path, dirName) + string(filepath.Separator)
	derr := os.MkdirAll(workingDir, 0777)
	if derr != nil {
		panic(derr)
	}
	mclWorkingDir := filepath.Join(workingDir, "mclOuput/")
	fmt.Println("MCL working dir: ", mclWorkingDir)
	derr = os.MkdirAll(mclWorkingDir, 0777)
	if derr != nil {
		panic(derr)
	}
	slog.Info("Working path becomes: " + workingDir)
	// change working dir cause mcl program dumps files into working dir
	os.Chdir(mclWorkingDir)
	abcFileName := fmt.Sprintf("%s/%s.abc", mclWorkingDir, exp.Name)
	convertTrajectoriesToAbcFormat(exp, abcFileName)
	tabFileName := fmt.Sprintf("%s/%s.tab", mclWorkingDir, exp.Name)
	mciFileName := fmt.Sprintf("%s/%s.mci", mclWorkingDir, exp.Name)
	mcxloadCmd := fmt.Sprintf("%smcxload", pathToMcl)
	cmd := exec.Command(mcxloadCmd, "-abc", abcFileName, "--stream-mirror", "-write-tab", tabFileName, "-o", mciFileName)
	var out bytes.Buffer
	var serr bytes.Buffer
	cmd.Stdout = &out
	cmd.Stderr = &serr
	err := cmd.Run()
	if err != nil {
		panic(err)
	}
	slog.Info("Output: " + out.String() + serr.String())
	// run the clusterings with different granularities
	for _, gran := range granularities {
		mcl_cmd := fmt.Sprintf("%smcl", pathToMcl)
		cmd := exec.Command(mcl_cmd, mciFileName, "-I", fmt.Sprintf("%f", float64(gran)/10.0))
		var out2 bytes.Buffer
		var serr2 bytes.Buffer
		cmd.Stdout = &out2
		cmd.Stderr = &serr2
		slog.Info("Output: " + out2.String() + serr2.String())
		err := cmd.Run()
		if err != nil {
			panic(err)
		}
	}
	// convert the clusterings to readable format
	clusterFileName := fmt.Sprintf("out.%s.mci", exp.Name)
	outFileName := fmt.Sprintf("%s.mci", exp.Name)
	mcxdumpCmd := fmt.Sprintf("%smcxdump", pathToMcl)
	for _, gran := range granularities {
		cmd := exec.Command(mcxdumpCmd, "-icl", fmt.Sprintf("%s.I%d", clusterFileName, gran), "-tabr", tabFileName, "-o", fmt.Sprintf("%s.I%d", outFileName, gran))
		slog.Info(fmt.Sprint(mcxdumpCmd, "-icl", fmt.Sprintf("%s.I%d", clusterFileName, gran), "-tabr", tabFileName, "-o", fmt.Sprintf("%s.I%d", outFileName, gran)))
		var out1 bytes.Buffer
		var serr1 bytes.Buffer
		cmd.Stdout = &out1
		cmd.Stderr = &serr1
		err := cmd.Run()
		slog.Info("Output: " + out1.String() + serr1.String())
		if err != nil {
			panic(err)
		}
	}
	// convert the clusterings generated by mcl tool to gml format
	for _, gran := range granularities {
		granDirName := fmt.Sprintf("%s-clusters-I%d/", exp.Name, gran)
		outPath := filepath.Join(workingDir, granDirName) + string(filepath.Separator)
		derr := os.MkdirAll(outPath, 0777)
		if derr != nil {
			panic(derr)
		}
		inputFileName := fmt.Sprintf("%s.I%d", outFileName, gran)
		convertToDirectTrajectoryClusterGraphs(exp, inputFileName, fmt.Sprintf("%s%s.trajectories.gml", outPath, inputFileName))
		convertToDirectTrajectoryClusterGraphsDot(exp, inputFileName, fmt.Sprintf("%s%s.trajectories.dot", outPath, inputFileName))
		convertToDirectTrajectoryClusterGraphsRR(exp, inputFileName, fmt.Sprintf("%s%s.trajectories.RR.gml", outPath, inputFileName))
		convertToDirectTrajectoryClusterGraphsRRDot(exp, inputFileName, fmt.Sprintf("%s%s.trajectories.RR.dot", outPath, inputFileName))
		trajectory.PrintClusteredTrajectoriesToFile(exp, fmt.Sprintf("%s%s.clustered.trajectories.tab", outPath, inputFileName))
		trajectory.PrintClustersToCSVFiles(exp, fmt.Sprintf("%s%s.clustered.patients.csv", outPath, inputFileName),
			fmt.Sprintf("%s%s.clustered.clusters.csv", outPath, inputFileName))
	}
}

// collectTrajectoriesFromClusterData looks up trajectories associated with a given list of trajectory ids and assigns
// each of these to a specific cluster id. It returns the list of trajectory objects.
func collectTrajectoriesFromClusterData(exp *trajectory.Experiment, ids []int, clusterID int) []*trajectory.Trajectory {
	ts := []*trajectory.Trajectory{}
	for _, id := range ids {
		// assign cluster label to trajectory
		exp.Trajectories[id].Cluster = clusterID
		ts = append(ts, exp.Trajectories[id])
	}
	return ts
}

// convertToDirectTrajectoryClusterGraphs produces a GML graph file for the clustered trajectories in an experiment. For
// this, it parses the cluster output from MCL, which is a file that lists for each cluster id a list of trajectory ids
// that are assigned to it. Then it looks up the concrete trajectory objects for each trajectory id. Finally, each
// cluster is written to the output file by writing all of the cluster's trajectories as part of a subgraph for that
// cluster.
func convertToDirectTrajectoryClusterGraphs(exp *trajectory.Experiment, input, output string) {
	file, err := os.Open(input)
	if err != nil {
		panic(err)
	}
	ofile, oerr := os.Create(output)
	if oerr != nil {
		panic(oerr)
	}
	defer func() {
		if err := file.Close(); err != nil {
			panic(err)
		}
		if oerr := ofile.Close(); oerr != nil {
			panic(oerr)
		}
	}()
	// trajectories to assign to clusters
	nofClusters := 0

	// parse file
	reader := csv.NewReader(file)
	reader.Comma = '\t'
	reader.FieldsPerRecord = -1
	reader.LazyQuotes = true
	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			panic(err)
		}
		// collect codes in the cluster
		var codes []int
		for _, rcode := range record {
			code, err := strconv.Atoi(rcode)
			if err != nil {
				panic(err)
			}
			codes = append(codes, code)
		}
		// print the trajectories in the cluster
		collected := collectTrajectoriesFromClusterData(exp, codes, nofClusters)
		nofClusters++
		// print this cluster
		// print header
		fmt.Fprintf(ofile, "graph [ \n directed 1 \n multigraph 1\n")
		nodePrinted := map[int]bool{}
		// print nodes
		for _, t := range collected {
			for _, node := range t.Diagnoses {
				if _, ok := nodePrinted[node]; !ok {
					fmt.Fprintf(ofile, fmt.Sprintf("node [ id %d\n label \"%s\"\n ]\n", node, exp.NameMap[node]))
					nodePrinted[node] = true
				}
			}
		}
		// print edges
		edgePrinted := make([][][]int, exp.NofDiagnosisCodes)
		for i, _ := range edgePrinted {
			edgePrinted[i] = make([][]int, exp.NofDiagnosisCodes)
		}
		for _, t := range collected {
			d1 := t.Diagnoses[0]
			for i := 1; i < len(t.Diagnoses); i++ {
				d2 := t.Diagnoses[i]
				n := t.PatientNumbers[i-1]
				printed := edgePrinted[d1][d2]
				if !utils.MemberInt(n, printed) {
					fmt.Fprintf(ofile, fmt.Sprintf("edge [\nsource %d\ntarget %d\nlabel %d\n]\n", d1, d2, n))
					if printed == nil {
						edgePrinted[d1][d2] = []int{n}
					} else {
						edgePrinted[d1][d2] = append(edgePrinted[d1][d2], n)
					}
				}
				d1 = d2
			}
		}
		fmt.Fprintf(ofile, "]\n")
	}
	slog.Info("For "+output,
		slog.Int("Clusters", nofClusters),
		slog.Int("Trajectories", len(exp.Trajectories)),
	)
}

// percentMalesFemales computes for a given list of patients the percentage of males and females wrt to the total number
// of males and females in the experiment.
func percentMalesFemales(exp *trajectory.Experiment, ps []*trajectory.Patient) (float64, float64) {
	m := 0
	f := 0
	for _, p := range ps {
		if p.Sex == trajectory.Male {
			m++
		} else {
			f++
		}
	}
	return (100.0 / float64(exp.MCtr)) * float64(m), (100.0 / float64(exp.FCtr)) * float64((f))
}

// getDiagnosisDate returns the concrete diagnosis date for a given pair of diagnosis ids.
func getDiagnosisDate(p *trajectory.Patient, d1, d2 int) trajectory.DiagnosisDate {
	d1idx := -1
	d2idx := -1
	for i, d := range p.Diagnoses {
		if d.DID == d1 {
			d1idx = i
			continue
		}
		if d.DID == d2 && d1idx != -1 {
			d2idx = i
		}
	}
	return p.Diagnoses[d2idx].Date
}

// percentEOI computes percent of patients that have their event of interest at the time of the transition of disease
// d1 -> d2
func percentEOI(exp *trajectory.Experiment, ps []*trajectory.Patient, d1, d2 int) float64 {
	eoictr := 0
	for _, p := range ps {
		d := getDiagnosisDate(p, d1, d2)
		if p.EOIDate != nil && trajectory.DiagnosisDateSmallerThan(*p.EOIDate, d) {
			eoictr++
		}
	}
	return (100.0 / float64(len(ps))) * float64(eoictr)
}

func transitionInformation(exp *trajectory.Experiment, t *trajectory.Trajectory, i, d1, d2 int) (string, string, string) {
	rr := strconv.FormatFloat(exp.DxDRR[d1][d2], 'f', 2, 64)
	m, f := percentMalesFemales(exp, t.Patients[i])
	mfratio := strconv.FormatFloat(m/f, 'f', 2, 64)
	eoi := strconv.FormatFloat(percentEOI(exp, t.Patients[i], d1, d2), 'f', 0, 64)
	return rr, mfratio, eoi
}

// convertToDirectTrajectoryClusterGraphsRR converts MCL cluster output - a file with for each cluster id a list of
// trajectory ids - to a GML output file that plots the trajectories as graphs. Each cluster is plotted as a separate
// subgraph, with diagnosis codes used as nodes and trajectory transitions used as edges. The edges are annotated with
// the relative risk score (RR) associated with the diagnosis pair that the edge represents.
func convertToDirectTrajectoryClusterGraphsRR(exp *trajectory.Experiment, input, output string) {
	file, err := os.Open(input)
	if err != nil {
		panic(err)
	}
	ofile, oerr := os.Create(output)
	if oerr != nil {
		panic(oerr)
	}
	defer func() {
		if err := file.Close(); err != nil {
			panic(err)
		}
		if oerr := ofile.Close(); oerr != nil {
			panic(oerr)
		}
	}()
	// trajectories to assign to clusters
	nofClusters := 0

	// parse file
	reader := csv.NewReader(file)
	reader.Comma = '\t'
	reader.FieldsPerRecord = -1
	reader.LazyQuotes = true
	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			panic(err)
		}
		// collect codes in the cluster
		var codes []int
		for _, rcode := range record {
			code, err := strconv.Atoi(rcode)
			if err != nil {
				panic(err)
			}
			codes = append(codes, code)
		}
		// print the trajectories in the cluster
		collected := collectTrajectoriesFromClusterData(exp, codes, nofClusters)
		nofClusters++
		// print this cluster
		// print header
		fmt.Fprintf(ofile,
			fmt.Sprintf("graph [ \n comment \"cluster %d\" \n directed 1 \n label \"cluster %d\" \n "+
				"multigraph 1\n", nofClusters-1, nofClusters-1))
		nodePrinted := map[int]bool{}
		// print nodes
		for _, t := range collected {
			for _, node := range t.Diagnoses {
				if _, ok := nodePrinted[node]; !ok {
					fmt.Fprintf(ofile, fmt.Sprintf("node [ id %d\n label \"%s\"\n ]\n", node, exp.NameMap[node]))
					nodePrinted[node] = true
				}
			}
		}
		// print edges
		edgePrinted := make([][]bool, exp.NofDiagnosisCodes)
		for i, _ := range edgePrinted {
			edgePrinted[i] = make([]bool, exp.NofDiagnosisCodes)
		}
		for _, t := range collected {
			d1 := t.Diagnoses[0]
			tctr := 0
			for i := 1; i < len(t.Diagnoses); i++ {
				d2 := t.Diagnoses[i]
				if !edgePrinted[d1][d2] {
					edgePrinted[d1][d2] = true
					RR := strconv.FormatFloat(exp.DxDRR[d1][d2], 'f', 2, 64)
					fmt.Fprintf(ofile, fmt.Sprintf("edge [\nsource %d\ntarget %d\nlabel %s\n]\n", d1, d2, RR))
					//rr, mfratio, eoi := transitionInformation(exp, t, tctr, d1, d2)
					//fmt.Fprintf(ofile, fmt.Sprintf("edge [\nsource %d\ntarget %d\nlabel \"RR:%s,M/F:%s,EOI:%s\"\n]\n", d1, d2, rr, mfratio, eoi))
				}
				d1 = d2
				tctr++
			}
		}
		fmt.Fprintf(ofile, "]\n")
	}
	slog.Info("For "+output,
		slog.Int("Clusters", nofClusters),
		slog.Int("Trajectories", len(exp.Trajectories)),
	)
}

// convertToDirectTrajectoryClusterGraphsRRDot converts MCL cluster output - a file with for each cluster id a list of
// trajectory ids - to a DOT output file that plots the trajectories as graphs using GraphViz. Each cluster is plotted
// as a separate subgraph, with diagnosis codes used as nodes and trajectory transitions used as edges. The edges are
// annotated with the relative risk score (RR) associated with the diagnosis pair that the edge represents.
func convertToDirectTrajectoryClusterGraphsRRDot(exp *trajectory.Experiment, input, output string) {
	file, err := os.Open(input)
	if err != nil {
		panic(err)
	}
	ofile, oerr := os.Create(output)
	if oerr != nil {
		panic(oerr)
	}
	defer func() {
		if err := file.Close(); err != nil {
			panic(err)
		}
		if oerr := ofile.Close(); oerr != nil {
			panic(oerr)
		}
	}()
	// trajectories to assign to clusters
	nofClusters := 0

	// parse file
	reader := csv.NewReader(file)
	reader.Comma = '\t'
	reader.FieldsPerRecord = -1
	reader.LazyQuotes = true

	// print dot directional-graph header, and some formatting
	fmt.Fprintf(ofile, "digraph {\n")
	fmt.Fprintf(ofile, "  node [fontsize=24 fillcolor=lightskyblue style=filled color=navy]")
	fmt.Fprintf(ofile, "  edge [fontcolor=red fontsize=30]")

	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			panic(err)
		}
		// collect codes in the cluster
		var codes []int
		for _, rcode := range record {
			code, err := strconv.Atoi(rcode)
			if err != nil {
				panic(err)
			}
			codes = append(codes, code)
		}
		// print the trajectories in the cluster
		collected := collectTrajectoriesFromClusterData(exp, codes, nofClusters)
		nofClusters++
		// print this cluster
		// print header
		fmt.Fprintf(ofile, fmt.Sprintf("  subgraph cluster_%d {\n", nofClusters-1))

		nodePrinted := map[int]bool{}
		// print nodes
		for _, t := range collected {
			for _, node := range t.Diagnoses {
				if _, ok := nodePrinted[node]; !ok {
					fmt.Fprintf(ofile, fmt.Sprintf("    c%d_%d [label=\"%s\"]\n", nofClusters-1, node, utils.WrapText(exp.NameMap[node], 25)))
					nodePrinted[node] = true
				}
			}
		}
		// print edges
		edgePrinted := make([][]bool, exp.NofDiagnosisCodes)
		for i, _ := range edgePrinted {
			edgePrinted[i] = make([]bool, exp.NofDiagnosisCodes)
		}
		for _, t := range collected {
			d1 := t.Diagnoses[0]
			tctr := 0
			for i := 1; i < len(t.Diagnoses); i++ {
				d2 := t.Diagnoses[i]
				if !edgePrinted[d1][d2] {
					edgePrinted[d1][d2] = true
					RR := strconv.FormatFloat(exp.DxDRR[d1][d2], 'f', 2, 64)
					label := RR + " [" + strconv.FormatFloat(exp.DxDRRPval[d1][d2], 'f', 2, 64) + "]"
					fmt.Fprintf(ofile,
						fmt.Sprintf("    c%d_%d -> c%d_%d [label=\"%s\" penwidth=%s weight=%s]\n", nofClusters-1, d1, nofClusters-1, d2, label, RR, RR))
				}
				d1 = d2
				tctr++
			}
		}
		fmt.Fprintf(ofile, "  }\n")
	}
	// close the DOT digraph
	fmt.Fprintf(ofile, "}\n")

	slog.Info("For "+output,
		slog.Int("Clusters", nofClusters),
		slog.Int("Trajectories", len(exp.Trajectories)),
	)
}

// convertToDirectTrajectoryClusterGraphsRRDot converts MCL cluster output - a file with for each cluster id a list of
// trajectory ids - to a DOT output file that plots the trajectories as graphs using GraphViz. Each cluster is plotted
// as a separate subgraph, with diagnosis codes used as nodes and trajectory transitions used as edges. The edges are
// annotated with the relative risk score (RR) associated with the diagnosis pair that the edge represents.
func convertToDirectTrajectoryClusterGraphsDot(exp *trajectory.Experiment, input, output string) {
	file, err := os.Open(input)
	if err != nil {
		panic(err)
	}
	ofile, oerr := os.Create(output)
	if oerr != nil {
		panic(oerr)
	}
	defer func() {
		if err := file.Close(); err != nil {
			panic(err)
		}
		if oerr := ofile.Close(); oerr != nil {
			panic(oerr)
		}
	}()
	// trajectories to assign to clusters
	nofClusters := 0

	// parse file
	reader := csv.NewReader(file)
	reader.Comma = '\t'
	reader.FieldsPerRecord = -1
	reader.LazyQuotes = true

	// print dot directional-graph header, and some formatting
	fmt.Fprintf(ofile, "digraph {\n")
	fmt.Fprintf(ofile, "  node [fontsize=24 fillcolor=lightskyblue style=filled color=navy]")
	fmt.Fprintf(ofile, "  edge [fontcolor=red fontsize=30]")

	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			panic(err)
		}
		// collect codes in the cluster
		var codes []int
		for _, rcode := range record {
			code, err := strconv.Atoi(rcode)
			if err != nil {
				panic(err)
			}
			codes = append(codes, code)
		}
		// print the trajectories in the cluster
		collected := collectTrajectoriesFromClusterData(exp, codes, nofClusters)
		nofClusters++
		// print this cluster
		// print header
		fmt.Fprintf(ofile, fmt.Sprintf("  subgraph cluster_%d {\n", nofClusters-1))

		nodePrinted := map[int]bool{}
		// print nodes
		for _, t := range collected {
			for _, node := range t.Diagnoses {
				if _, ok := nodePrinted[node]; !ok {
					fmt.Fprintf(ofile, fmt.Sprintf("    c%d_%d [label=\"%s\"]\n", nofClusters-1, node, utils.WrapText(exp.NameMap[node], 25)))
					nodePrinted[node] = true
				}
			}
		}
		// print edges
		edgePrinted := make([][]bool, exp.NofDiagnosisCodes)
		for i, _ := range edgePrinted {
			edgePrinted[i] = make([]bool, exp.NofDiagnosisCodes)
		}
		for _, t := range collected {
			d1 := t.Diagnoses[0]
			tctr := 0
			for i := 1; i < len(t.Diagnoses); i++ {
				d2 := t.Diagnoses[i]
				if !edgePrinted[d1][d2] {
					edgePrinted[d1][d2] = true
					v := float64(t.PatientNumbers[i-1]) / 10.0
					n := strconv.FormatInt(int64(t.PatientNumbers[i-1]), 10)
					ns := strconv.FormatFloat(v, 'f', 2, 64)
					label := n + " [" + strconv.FormatFloat(exp.DxDRRPval[d1][d2], 'f', 2, 64) + "]"
					fmt.Fprintf(ofile,
						fmt.Sprintf("    c%d_%d -> c%d_%d [label=\"%s\" penwidth=%s weight=%s]\n", nofClusters-1, d1, nofClusters-1, d2, label, ns, ns))
				}
				d1 = d2
				tctr++
			}
		}
		fmt.Fprintf(ofile, "  }\n")
	}
	// close the DOT digraph
	fmt.Fprintf(ofile, "}\n")

	slog.Info("For "+output,
		slog.Int("Clusters", nofClusters),
		slog.Int("Trajectories", len(exp.Trajectories)),
	)
}

// clusterGraph is a graph representation of a list of trajectory clusters.
type clusterGraph struct {
	nodes             []int         //list of nodes in the graph. A node is a diagnosis code that occurs in one or more trajectories
	edges             map[int][]int //adjacency list of edges
	clusterMembership map[int][]int //per node/diagnosis code the clusters it belongs to
	weight            int
	degrees           map[int]int
	similarities      map[int]map[int]float64
}

// createClusterGraph creates a new clusterGraph object
func createClusterGraph() clusterGraph {
	//clusterMembership and graph can only be initialised when the nodes
	return clusterGraph{nodes: []int{}, edges: map[int][]int{},
		clusterMembership: map[int][]int{},
		degrees:           map[int]int{},
		similarities:      map[int]map[int]float64{}}
}

// addClusterMembership registers that a diagnosis code occurs in a specific cluster.
func (graph *clusterGraph) addClusterMembership(did int, cid int) {
	entries, ok := graph.clusterMembership[did]
	if !ok {
		entries = []int{cid}
		graph.clusterMembership[did] = entries
		return
	}
	if !(utils.MemberInt(cid, entries)) { //add membership only once
		graph.clusterMembership[did] = append(entries, cid)
	}
}

// addEdge adds a diagnosis pair to the graph as a unique edge.
func (graph *clusterGraph) addEdge(d1 int, d2 int) {
	entries, ok := graph.edges[d1]
	if !ok {
		entries = []int{d2}
		graph.edges[d1] = entries
		return
	}
	if !(utils.MemberInt(d2, entries)) { //add edge only once
		graph.edges[d1] = append(entries, d2)
	}
}

// addNode adds a diagnosis code to the graph as a unique node.
func (graph *clusterGraph) addNode(d int) {
	ok := utils.MemberInt(d, graph.nodes) //add node only once
	if !ok {
		graph.nodes = append(graph.nodes, d)
		graph.degrees[d] = 0
	}
}

// addTrajectories adds all nodes and edges from a list of trajectories to the graph
func (graph *clusterGraph) addTrajectories(trajectories []*trajectory.Trajectory, cid int) {
	for _, t := range trajectories {
		d1 := t.Diagnoses[0]
		graph.addNode(d1)
		graph.addClusterMembership(d1, cid)
		for _, d2 := range t.Diagnoses[1:] {
			graph.addNode(d2)
			graph.addClusterMembership(d2, cid)
			graph.addEdge(d1, d2)
			d1 = d2
		}
	}
}

// calculateWeightsAndDegrees calculates:
// - the graph weight = the total number of edges
// - the node degrees = the number of edges per node
func (graph *clusterGraph) calculateWeightsAndDegrees() {
	//calculate weight
	for d, ds := range graph.edges {
		graph.weight += len(ds)
		_, ok := graph.degrees[d]
		if !ok {
			graph.degrees[d] = len(ds)
		} else {
			graph.degrees[d] += len(ds)
		}
	}
}

func (graph *clusterGraph) clustersShared(d1, d2 int) int {
	c1s := graph.clusterMembership[d1]
	c2s := graph.clusterMembership[d2]
	ctr := 0
	for _, c1 := range c1s {
		if utils.MemberInt(c1, c2s) {
			ctr++
		}
	}
	return ctr
}

func (graph *clusterGraph) nClusters(d int) int {
	v, ok := graph.clusterMembership[d]
	if !ok {
		panic("Disease not part of cluster. Should not happen.")
	}
	return len(v)
}

func (graph *clusterGraph) calculateEdgeSimilarities() {
	for d1, ds := range graph.edges {
		if _, ok := graph.similarities[1]; !ok {
			graph.similarities[d1] = map[int]float64{} //init similarity scores
		}
		for _, d2 := range ds {
			sharedN := graph.clustersShared(d1, d2)
			d1N := graph.nClusters(d1)
			d2N := graph.nClusters(d2)
			s := float64(sharedN) / float64(d1N*d2N)
			graph.similarities[d1][d2] = s
		}
	}
}

func (graph *clusterGraph) print() {
	fmt.Println("Graph: ")
	fmt.Println("Nodes: ")
	fmt.Println(graph.nodes)
	fmt.Println("Degrees: ")
	fmt.Println(graph.degrees)
	fmt.Println("Edges: ")
	fmt.Println(graph.edges)
	fmt.Println("Weight: ")
	fmt.Println(graph.weight)
	fmt.Println("Similarties: ")
	fmt.Println(graph.similarities)
	fmt.Println("Cluster membership: ")
	fmt.Println(graph.clusterMembership)
}

func (graph *clusterGraph) ExtendedModularityMetric() float64 {
	graph.calculateWeightsAndDegrees()
	graph.calculateEdgeSimilarities()
	m2 := float64(2 * graph.weight)
	contributions := 0.0
	for d1, ds := range graph.edges {
		for _, d2 := range ds {
			w := 1.0
			degree1 := graph.degrees[d1]
			degree2 := graph.degrees[d2]
			e := float64(degree1*degree2) / m2
			contributions += (w - e) * graph.similarities[d1][d2]
		}
	}
	return contributions / m2
}

// parseMCLGraph takes as input a file with a mcl clustering. The mcl file contains per cluster ID a list of trajectory
// IDs that belong to that cluster. Create an clusterGraph object with the following information that can be used to
// calculate an extended modularity metric of a graph:
// the list of diagnosis codes in the mcl graph
// the list of edges in the mcl graph
// per diagnosis code the list of clusters it occurs in
func ParseMCLGraph(exp *trajectory.Experiment, graph string) clusterGraph {

	file, err := os.Open(graph)
	if err != nil {
		panic(err)
	}
	defer func() {
		if err := file.Close(); err != nil {
			panic(err)
		}
	}()

	// trajectories to assign to clusters
	nofClusters := 0

	// parse file
	reader := csv.NewReader(file)
	reader.Comma = '\t'
	reader.FieldsPerRecord = -1
	reader.LazyQuotes = true

	mclGraph := createClusterGraph()

	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			panic(err)
		}
		// collect the trajectory codes in the cluster
		var tids []int
		for _, c := range record {
			tid, err := strconv.Atoi(c)
			if err != nil {
				panic(err)
			}
			tids = append(tids, tid)
		}
		// print the trajectories in the cluster
		trajectories := collectTrajectoriesFromClusterData(exp, tids, nofClusters)
		//add the trajectory info to the graph
		mclGraph.addTrajectories(trajectories, nofClusters)
		nofClusters++
	}
	return mclGraph
}
