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
	"os"
	"os/exec"
	"path/filepath"
	"ptra/trajectory"
	"ptra/utils"
	"strconv"
)

// Clustering as in Brunak paper
// Cluster all diagnostic codes in the computed trajectories. This results in clusters of diagnostic codes. Clusters are
// visualised by plotting all trajectories that have all their diagnostic codes in that cluster as directed graphs.
// The trajectories are seen as sets of diagnostic codes that likely have overlap. The diagnostic codes can be modeled
// together as a graph. The connection between diagnostic codes represents a diagnosis pair. The strength between
// diagnostic codes is measured with the jaccard similarity coefficient. This is computed as:
// total nr of trajectories (sets) a pair A->B occurs in / total trajectories A occurs in + total trajectories B
// occurs in - total nr of trajectories A->B occurs in. The jaccard index is a normalized value [0, 1] for all pairs.
// Doing the clustering this way preserves the importance of the directionality of the pairs.

// computeTotalOccurencesPairs computes for every pair A->B detected:
// the total number of trajectories A belongs to
// the total number of trajectories B belongs to
// the total number of trajectories A->B belongs to
func computeTotalOccurencesPairs(exp *trajectory.Experiment) ([]int, [][]int) {
	diagnosisCounts := make([]int, exp.NofDiagnosisCodes)
	pairCounts := make([][]int, exp.NofDiagnosisCodes)
	for i, _ := range pairCounts {
		pairCounts[i] = make([]int, exp.NofDiagnosisCodes)
	}
	for _, t := range exp.Trajectories {
		d1 := t.Diagnoses[0]
		diagnosisCounts[d1]++
		for j := 1; j < len(t.Diagnoses); j++ {
			d2 := t.Diagnoses[j]
			diagnosisCounts[d2]++
			pairCounts[d1][d2]++
			d1 = d2
		}
	}
	return diagnosisCounts, pairCounts
}

// computeJaccardIndexForPairs computes for each diagnosis pair A->B the jaccard similarity coefficient.
func computeJaccardIndexForPairs(exp *trajectory.Experiment) [][]float64 {
	//create and initialise index. Jaccard index -1.0 means pair does not exist.
	index := make([][]float64, exp.NofDiagnosisCodes)
	for i, _ := range index {
		row := make([]float64, exp.NofDiagnosisCodes)
		for j, _ := range row {
			row[j] = -1.0
		}
		index[i] = row
	}
	diagnosisCounts, pairCounts := computeTotalOccurencesPairs(exp)
	for _, pair := range exp.Pairs {
		pairTotal := float64(pairCounts[pair.First][pair.Second])
		firstTotal := float64(diagnosisCounts[pair.First])
		secondTotal := float64(diagnosisCounts[pair.Second])
		jaccardCoeff := pairTotal / (firstTotal + secondTotal - pairTotal)
		index[pair.First][pair.Second] = jaccardCoeff
	}
	return index
}

func convertTrajectoryPairsToAbcFormat(exp *trajectory.Experiment, name string) {
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

	// compute the jacard index for the pairs
	jaccardIndex := computeJaccardIndexForPairs(exp)
	// plot pairs as part of the same graph
	for d1, d2s := range jaccardIndex {
		for d2, coeff := range d2s {
			if coeff >= 0 {
				fmt.Fprintf(file, "%d\t%d\t%f\n", d1, d2, coeff)
			}
		}
	}
}

func ClusterTrajectories(exp *trajectory.Experiment, granularities []int, path, pathToMcl string) {
	fmt.Println("Clustering trajectories with MCL")
	// convert trajectories to abc format for the mcl tool
	dirName := fmt.Sprintf("%s-clusters/", exp.Name)
	workingDir := filepath.Join(path, dirName) + string(filepath.Separator)
	fmt.Println("Working path becomes: ", workingDir)
	derr := os.MkdirAll(workingDir, 0777)
	if derr != nil {
		panic(derr)
	}
	// change working dir cause mcl program dumps files into working dir
	os.Chdir(workingDir)
	abcFileName := fmt.Sprintf("%s%s.abc", workingDir, exp.Name)
	convertTrajectoryPairsToAbcFormat(exp, abcFileName)
	tabFileName := fmt.Sprintf("%s%s.tab", workingDir, exp.Name)
	mciFileName := fmt.Sprintf("%s%s.mci", workingDir, exp.Name)
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
	fmt.Println("Output: ", out.String(), serr.String())
	// run the clusterings with different granularities
	for _, gran := range granularities {
		mcl_cmd := fmt.Sprintf("%smcl", pathToMcl)
		cmd := exec.Command(mcl_cmd, mciFileName, "-I", fmt.Sprintf("%f", float64(gran)/10.0))
		var out2 bytes.Buffer
		var serr2 bytes.Buffer
		cmd.Stdout = &out2
		cmd.Stderr = &serr2
		fmt.Println("Output: ", out2.String(), serr2.String())
		err := cmd.Run()
		if err != nil {
			panic(err)
		}
	}
	// convert the clusterings to readable format
	clusterFileName := fmt.Sprintf("out.%s.mci", exp.Name)
	outFileName := fmt.Sprintf("dump.%s.mci", exp.Name)
	mcxdumpCmd := fmt.Sprintf("%smcxdump", pathToMcl)
	for _, gran := range granularities {
		cmd := exec.Command(mcxdumpCmd, "-icl", fmt.Sprintf("%s.I%d", clusterFileName, gran), "-tabr", tabFileName, "-o", fmt.Sprintf("%s.I%d", outFileName, gran))
		fmt.Println(mcxdumpCmd, "-icl", fmt.Sprintf("%s.I%d", clusterFileName, gran), "-tabr", tabFileName, "-o", fmt.Sprintf("%s.I%d", outFileName, gran))
		var out1 bytes.Buffer
		var serr1 bytes.Buffer
		cmd.Stdout = &out1
		cmd.Stderr = &serr1
		err := cmd.Run()
		fmt.Println("Output: ", out1.String(), serr1.String())
		if err != nil {
			panic(err)
		}
	}
	// convert the clusterings generated by mcl tool to gml format
	for _, gran := range granularities {
		dumpFileName := fmt.Sprintf("%s.I%d", outFileName, gran)
		convertToTrajectoryClusterGraphs(exp, dumpFileName, fmt.Sprintf("%s.trajectories.gml", dumpFileName))
		convertToDiagnosisGraphs(exp, dumpFileName, fmt.Sprintf("%s.gml", dumpFileName))
	}
}

// collectTrajectoriesInCluster collects all trajectories that have all diagnosis codes in the cluster. Allow n missing
// diagnoses. (Brunak paper allows 1 miss)
func collectTrajectoriesInCluster(trajectories []*trajectory.Trajectory, cluster []int, n int) ([]*trajectory.Trajectory, []*trajectory.Trajectory) {
	collected := []*trajectory.Trajectory{}
	uncollected := []*trajectory.Trajectory{}
	for _, t := range trajectories {
		misses := 0
		pass := true
		for _, d := range t.Diagnoses {
			if !utils.MemberInt(d, cluster) {
				misses++
				if misses > n {
					pass = false
					continue
				}
			}
		}
		if pass {
			collected = append(collected, t)
		} else {
			uncollected = append(uncollected, t)
		}
	}
	return collected, uncollected
}

func convertToDiagnosisGraphs(exp *trajectory.Experiment, input, output string) {
	in, err := os.Open(input)
	if err != nil {
		panic(err)
	}
	out, err := os.Create(output)
	if err != nil {
		panic(err)
	}
	defer func() {
		if err := in.Close(); err != nil {
			panic(err)
		}
		if oerr := out.Close(); oerr != nil {
			panic(oerr)
		}
	}()
	// parse file
	reader := csv.NewReader(in)
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
		// print nodes
		fmt.Fprintf(out, "graph [ \n directed 1 \n multigraph 1\n")
		// print nodes
		for _, code := range codes {
			fmt.Fprintf(out, fmt.Sprintf("node [ id %d\n label \"%s\"\n ]\n", code, exp.NameMap[code]))
		}
		// print edges, i.e. for every node combo, print an edge if there exists a pair
		existingPairs := map[int]map[int]bool{}
		for _, p := range exp.Pairs {
			if ff, ok := existingPairs[p.First]; !ok {
				f := map[int]bool{}
				f[p.Second] = true
				existingPairs[p.First] = f
			} else {
				ff[p.Second] = true
			}
		}
		for _, d1 := range codes {
			for _, d2 := range codes {
				if f, ok := existingPairs[d1]; ok {
					if _, ok2 := f[d2]; ok2 {
						fmt.Fprintf(out, fmt.Sprintf("edge [\nsource %d\ntarget %d\n]\n", d1, d2))
					}
				}
			}
		}
		fmt.Fprintf(out, "]\n")
	}
}

// concertToTrajectoryClusters converts a MCI file to a trajectory cluster. The MCI file contains per line a cluster. The
// line lists all nodes/diagnosis codes that belong to to that cluster.
// We collect the trajectories that are fully contained in those clusters and plot them as a directed graph.
func convertToTrajectoryClusterGraphs(exp *trajectory.Experiment, input, output string) {
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
	trajectories := exp.Trajectories
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
		// collect the trajectories in the cluster
		collected, uncollected := collectTrajectoriesInCluster(trajectories, codes, 1)
		trajectories = uncollected
		if len(collected) > 0 {
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
	}
	// print the unclustered trajectories as separate clusters
	for _, t := range trajectories {
		fmt.Fprintf(ofile, "graph [ \n directed 1 \n multigraph 1\n")
		// print nodes
		for _, d := range t.Diagnoses {
			fmt.Fprintf(ofile, fmt.Sprintf("node [ id %d\n label \"%s\"\n ]\n", d, exp.NameMap[d]))
		}
		// print edges
		d1 := t.Diagnoses[0]
		for i := 1; i < len(t.Diagnoses); i++ {
			d2 := t.Diagnoses[i]
			n := t.PatientNumbers[i-1]
			fmt.Fprintf(ofile, fmt.Sprintf("edge [\nsource %d\ntarget %d\nlabel %d\n]\n", d1, d2, n))
			d1 = d2
		}
		fmt.Fprintf(ofile, "]\n")
	}
	fmt.Println("For ", output)
	fmt.Println("Collected ", nofClusters, " clusters and ", len(trajectories), " not clustered trajectories.")
	fmt.Println("Clustered ", len(exp.Trajectories)-len(trajectories), " out of ", len(exp.Trajectories), " trajectories.")
}
