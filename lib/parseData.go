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
	"encoding/csv"
	"encoding/json"
	"encoding/xml"
	"fmt"
	"github.com/imec-int/ptra/lib/utils"
	"io"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
	"time"
)

//Package ptra implements a patient trajectory analysis tool.
//The ptra program has 3 data inputs:
//A file mapping diagnosis ID (DID) -> medical Name.
//A file with patient info, associating a patient ID (PID) with YOB, sex, etc.
//A file with diagnosis info, mapping PID -> DID, date.

//TriNetX Database stores diagnoses data using a mix of ICD10 and ICD9 codes.
//We have an additional file that maps ICD9 IDs -> ICD10 IDs.
//We can download the ICD10 ID -> medical Name from https://www.cms.gov/medicare/icd-10/2022-icd-10-cm as an xml file.
//TriNetX stores patient info as a csv file, as well as the diagnoses info.
//For the medical Name mapping, we can also use the ICD10 -> CCSR Mapping which maps ICD10 onto 530 Categories with medical meaning. This
//mapping can be downloaded from https://www.hcup-us.ahrq.gov/toolssoftware/ccsr/dxccsr.jsp#download as a CSV file. This
//mapping performs a mapping ICD10 -> CCSR Categories -> medical Name.

// Parsing ICD10 hierarchy from xml
// Structs for unmarshalling ICD10 xml data
// Structure of an ICD10 code: ABC.XYZ(D): up to 7 characters
// The xml file structures these codes using:
// <chapter> <section> <desc> </desc> </section> </chapter>
// with chapter the first Level of the diagnosis code (A), section the second Level (D) and desc the rest (C.XYZ(D))

// icd10Diag captures the lowest levels of the ICD10 code (diag)
type icd10Diag struct {
	Name      string      `xml:"name"` //Unique diagnosis ID (DID) in ICD10 encoding
	Desc      string      `xml:"desc"` //A medical Name/description for a DID
	Diagnoses []icd10Diag `xml:"diag"` //A diagnosis can be split into more detailed diagnoses descriptors.
}

// icd10Section captures the second Level of the ICD10 code (section)
type icd10Section struct {
	Desc      string      `xml:"desc"`    //A medical Name/description for a DID.
	Id        string      `xml:"id,attr"` //Unique diagnosis ID (DID) in ICD10 encoding
	Diagnoses []icd10Diag `xml:"diag"`    //A diagnosis can be split into more detailed diagnoses descriptors.
}

// icd10Chapter captures the first (highest) Level of the ICD10 code (chapter)
type icd10Chapter struct {
	XmlName  xml.Name       `xml:"chapter"`
	Desc     string         `xml:"desc"`
	Sections []icd10Section `xml:"section"`
}

// icd10Hierarchy contains the full xml table with the ICD10 code hierarchy.
type icd10Hierarchy struct {
	XmlName  xml.Name       `xml:"ICD10CM.tabular"`
	Chapters []icd10Chapter `xml:"chapter"`
}

// parseIcd10HierarchyFromXML parses the xml file with the ICD10 hierarchy into an icd10Hierarchy object.
func parseIcd10HierarchyFromXml(file string) icd10Hierarchy {
	fmt.Println("Parsing ICD10 code hierarchy from XML file: ", file)
	//open file
	xmlFile, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer xmlFile.Close()
	xmlFileBytes, _ := io.ReadAll(xmlFile)
	//unmarshall
	icd10Hierarchy := icd10Hierarchy{}
	xml.Unmarshal(xmlFileBytes, &icd10Hierarchy)
	return icd10Hierarchy
}

// printIcd10Hierarchy prints an ICD10 hierarchy parsed from an XML file.
func printIcd10Hierarchy(hierarchy icd10Hierarchy) {
	fmt.Println("Printing ICD10 code hierarchy.")
	// count # DID per Level
	ctr1, ctr2, ctr3, ctr4, ctr5, ctr6, ctr7 := 0, 0, 0, 0, 0, 0, 0
	for _, chap := range hierarchy.Chapters {
		// Level 1
		ctr1++
		fmt.Println("Chapter: ", chap.Desc)
		for _, section := range chap.Sections {
			// Level 2
			ctr2++
			fmt.Println("Section: ", section.Desc)
			for _, diag := range section.Diagnoses {
				// Level 3
				ctr3++
				fmt.Println(diag.Name, " : ", diag.Desc)
				if len(diag.Diagnoses) == 0 {
					continue
				}
				for _, diag := range diag.Diagnoses {
					// Level 4
					ctr4++
					fmt.Println(diag.Name, " : ", diag.Desc)
					if len(diag.Diagnoses) == 0 {
						continue
					}
					for _, diag := range diag.Diagnoses {
						// Level 5
						ctr5++
						fmt.Println(diag.Name, " : ", diag.Desc)
						if len(diag.Diagnoses) == 0 {
							continue
						}
						for _, diag := range diag.Diagnoses {
							// Level 6
							ctr6++
							fmt.Println(diag.Name, " : ", diag.Desc)
							if len(diag.Diagnoses) == 0 {
								continue
							}
							// Level 7
							ctr7++
							fmt.Println(diag.Name, " : ", diag.Desc)
						}
					}
				}
			}
		}
	}
	fmt.Println("#ICD10 codes/descriptors per Level: ")
	fmt.Println("Lvl 0: ", ctr1, " Lvl 1: ", ctr2, " Lvl 2: ",
		ctr3, " Lvl 3: ", ctr4, " Lvl 4: ", ctr5, " Lvl 5: ", ctr6, " Lvl 6: ", ctr7)
}

//The ptra program needs a names map that maps DID -> medical Name. The following code extracts a Name map from an ICD10
//hierarchy and a given Level.

// Icd10Entry is a struct for containing a medical Name + Level + the Categories of a DID in ICD10 encoding.
type Icd10Entry struct {
	Name       string    // medical Name for a DID in ICD10 encoding
	Categories [6]string // the names of the ICD10 encoding higher and lower in the hierarchy.
	Level      int       // the ICD10 hierarchy Level of this Name.
}

type icd10Table struct {
	namesMap map[string]Icd10Entry //maps ICD10 DID to a medical Name, Level, and Categories to which it belongs.
}

// selectParentCategory returns the ICD10 Name of the parent category to which a DID belongs.
func selectParentCategory(name Icd10Entry) string {
	if name.Level > 0 {
		return name.Categories[name.Level-1]
	}
	return "None"
}

func printIcd10NameMap(table map[string]Icd10Entry) {
	fmt.Println("ICD10 Name map: ")
	for id, name := range table {
		fmt.Println(id, " : ", name)
	}
}

// initializeIcd10NameMap initializes a Name map for ICD10 DID -> medical Name, Level, and Categories it belongs to.
func initializeIcd10NameMap(file string) map[string]Icd10Entry {
	icd10Map := map[string]Icd10Entry{} // maps ICD10 DID to a medical Name, Level, and Categories to which it belongs.
	icd10Hierarchy := parseIcd10HierarchyFromXml(file)
	for _, chap := range icd10Hierarchy.Chapters {
		category0 := chap.Desc
		for _, section := range chap.Sections {
			category1 := section.Desc
			// manually unrolled loop since we know hierarchy is max 7 levels, otherwise recursive code
			for _, diag := range section.Diagnoses {
				if len(diag.Diagnoses) == 0 {
					icd10Name := Icd10Entry{
						Name:       diag.Desc,
						Categories: [6]string{category0, category1, "NONE", "NONE", "NONE", "NONE"},
						Level:      2,
					}
					icd10Map[diag.Name] = icd10Name
					continue
				}
				category2 := diag.Desc
				for _, diag := range diag.Diagnoses {
					if len(diag.Diagnoses) == 0 {
						icd10Name := Icd10Entry{
							Name:       diag.Desc,
							Categories: [6]string{category0, category1, category2, "NONE", "NONE", "NONE"},
							Level:      3,
						}
						icd10Map[diag.Name] = icd10Name
						continue
					}
					category3 := diag.Desc
					for _, diag := range diag.Diagnoses {
						if len(diag.Diagnoses) == 0 {
							ICD10Name := Icd10Entry{
								Name:       diag.Desc,
								Categories: [6]string{category0, category1, category2, category3, "NONE", "NONE"},
								Level:      4,
							}
							icd10Map[diag.Name] = ICD10Name
							continue
						}
						category4 := diag.Desc
						for _, diag := range diag.Diagnoses {
							if len(diag.Diagnoses) == 0 {
								ICD10Name := Icd10Entry{
									Name:       diag.Desc,
									Categories: [6]string{category0, category1, category2, category3, category4, "NONE"},
									Level:      5,
								}
								icd10Map[diag.Name] = ICD10Name
								continue
							}
							category5 := diag.Desc
							for _, diag := range diag.Diagnoses {
								ICD10Name := Icd10Entry{
									Name:       diag.Desc,
									Categories: [6]string{category0, category1, category2, category3, category4, category5},
									Level:      6,
								}
								icd10Map[diag.Name] = ICD10Name
							}
						}
					}
				}
			}
		}
	}
	return icd10Map
}

// getIcd10DescToExcludeFromAnalysis returns a map that lists ICD10 Categories to be excluded from analysis by mapping
// the ICD10 category description (string) onto a boolean.
func getIcd10DescToExcludeFromAnalysis() map[string]bool {
	exclude := map[string]bool{}
	exclude["Pregnancy, childbirth and the puerperal (O00-O9A)"] = true
	exclude["Certain conditions originating in the perinatal period (P00-P96)"] = true
	exclude["Symptoms, signs and abnormal clinical and laboratory findings, not elsewhere classified (R00-R99)"] = true
	exclude["Injury, poisoning and certain other consequences of external causes (S00-T88)"] = true
	exclude["External causes of morbidity (V00-Y99)"] = true
	exclude["Factors influencing health status and contact with health services (Z00-Z99)"] = true
	return exclude
}

// getIcd10CodesToExcludeFromAnalysis returns the first letters of ICD10 codes to exclude from analysis.
func getIcd10CodesToExcludeFromAnalysis() map[string]bool {
	exclude := map[string]bool{}
	exclude["O"] = true
	exclude["P"] = true
	exclude["R"] = true
	exclude["S"] = true
	exclude["T"] = true
	exclude["V"] = true
	exclude["X"] = true
	exclude["Y"] = true
	exclude["Z"] = true
	return exclude
}

// getNonICD10CodesToAddToAnalysis returns a set of mockup ICD10 codes to be able to introduce non ICD codes to be
// included for analysis. It returns a map from mockup ICD10 code (string) to description string. It introduces "C98" for
// "Radical custectomy (bladder cancer)", "C99" for "MVAC Chemotherapy (bladder cancer)", and "C100" for "Intravesical
// therapy (bladder cancer)".
func getNonICD10CodesToAddToAnalysis() map[string]string {
	return map[string]string{
		"C98":  "Radical cystectomy (bladder cancer)",
		"C99":  "MVAC Chemotherapy (bladder cancer)",
		"C100": "Intravesical therapy (bladder cancer)",
	}
}

// initializeIcd10AnalysisIDMap creates a map ICD10 DID -> analysis DID and a map analysis ID -> medical Name. This is
// useful to remap diagnosis codes used in the input to a higher Level in the ICD10 hierarchy. E.g "typhoid fever" and
// "cholera" are both "infectious intestinal diseases", so they could both be identified as such during the analysis.
// This can be interesting to obtain more global patient trajectories/clusters.
func initializeIcd10AnalysisMaps(icd10Map map[string]Icd10Entry, level int) (map[string]int, map[int]Icd10Entry, int) {
	analysisIdMap := map[string]int{}                     // maps icd 10 code to analysis ID
	analysisIcd10Map := map[int]Icd10Entry{}              // maps analysis ID to an Icd10Entry
	nameToAnalysisIdMap := map[string]int{}               // maps medical Name to analysis ID
	ctr := 0                                              //serves as analysis ID generator
	icd10ToExclude := getIcd10DescToExcludeFromAnalysis() // a list of Level 0 Categories to exclude from analysis
	for icd10Code, icd10Entry := range icd10Map {
		if _, ok := icd10ToExclude[icd10Entry.Categories[0]]; ok {
			// code to exclude from analysis
			continue
		}
		var name string
		if level >= icd10Entry.Level {
			name = icd10Entry.Name
		} else {
			name = icd10Entry.Categories[level]
		}
		// may already have seen Name, because of Level
		newID, ok := nameToAnalysisIdMap[name]
		if !ok {
			newID = ctr
			ctr++
			analysisIcd10Map[newID] = icd10Entry
			nameToAnalysisIdMap[name] = newID
		}
		analysisIdMap[icd10Code] = newID
	}
	extra := getNonICD10CodesToAddToAnalysis()
	for code, name := range extra {
		analysisIcd10Map[ctr] = Icd10Entry{Name: name}
		nameToAnalysisIdMap[name] = ctr
		analysisIdMap[code] = ctr
		ctr++
	}
	fmt.Println("Mapped ", len(icd10Map), " ICD10 codes to ", ctr, " analysis IDs of Level ", level)
	return analysisIdMap, analysisIcd10Map, ctr
}

// ccsrCategory is a struct for containing CCSR categories, encoding medically meaningful names for a DID in ICD10
// encoding.
type ccsrCategory struct {
	name       string            //default CCSR category/medical name
	id         string            //CCSR ID for default category
	categories map[string]string //Up to 6 different CCSR categories an ICD10 code is mapped to
}

type icd10ToCCSRTable map[string]ccsrCategory //maps ICD10 DID to its CCSR Categories

// ccsrIcd10ToProperIcd10 transforms the ICD10 code from a ccsr file into a proper ICD10 code. The ICD10 codes in the
// ccsr file are stored without the ".", so this needs to be added to be able to compare to any other data that uses
// ICD10 codes. Also removes superfluous quotes stored in the ccsr file.
func ccsrIcd10ToProperIcd10(code string) string {
	return code[1:4] + "." + code[4:len(code)-1]
}

// initializeIcd10NameMapFromCCSR initializes a Name map for ICD10 DID -> CCSR Categories (medical names)
func initializeIcd10ToCCSRMap(file string) map[string]ccsrCategory {
	//map to collect data
	icd10ToCCSRTable := map[string]ccsrCategory{}
	//open file
	csvFile, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer func() {
		if err := csvFile.Close(); err != nil {
			panic(err)
		}
	}()
	//parse file
	reader := csv.NewReader(csvFile)
	//the header is 'ICD-10-CM CODE','ICD-10-CM CODE DESCRIPTION','Default CCSR CATEGORY IP','
	//Default CCSR CATEGORY DESCRIPTION IP','Default CCSR CATEGORY OP','Default CCSR CATEGORY DESCRIPTION OP','
	//CCSR CATEGORY 1','CCSR CATEGORY 1 DESCRIPTION','CCSR CATEGORY 2','CCSR CATEGORY 2 DESCRIPTION',
	//'CCSR CATEGORY 3','CCSR CATEGORY 3 DESCRIPTION','CCSR CATEGORY 4','CCSR CATEGORY 4 DESCRIPTION',
	//'CCSR CATEGORY 5','CCSR CATEGORY 5 DESCRIPTION','CCSR CATEGORY 6','CCSR CATEGORY 6 DESCRIPTION'
	// skip header
	reader.Read()
	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			panic(err)
		}
		//create CSSR category, set default category
		category := ccsrCategory{name: record[2], id: record[3], categories: map[string]string{}}
		//fill in unique CSSR alternative Categories, up to 6 possible
		for i := 6; i <= 17; i = i + 2 {
			if i >= len(record) {
				panic("Invalid DiagnosisInfo file")
			}
			catID := record[i]
			catName := record[i+1]
			if catName == "" || catID == "' '" {
				continue
			}
			if _, ok := category.categories[catID]; !ok {
				category.categories[catID] = catName
			}
		}
		//add category to result
		icd10Code := ccsrIcd10ToProperIcd10(record[0])
		icd10ToCCSRTable[icd10Code] = category
	}
	return icd10ToCCSRTable
}

// printIcd10ToCSSRTable is a simple function to print the map from iCD10 code to ccsr category. Useful for debugging.
func printIcd10ToCCSRTable(tab map[string]ccsrCategory) {
	fmt.Println("ICD10 to CCSR table")
	ctr := 0
	for icd10Code, ccsr := range tab {
		ctr++
		fmt.Println(icd10Code, " : ", ccsr.categories)
		if ctr >= 1000 {
			return
		}
	}
}

// initializeIcd10AnalysisMapsCCSR creates a map ICD10 DID -> [analysis DID] and a map analysis ID -> medical Name,
// starting from a CCSR mapping, which maps ICD10 codes onto medical meaningful Categories.
// Each icd10 code can be mapped to multiple ccsr Categories, and therefore to multiple analysis IDs.
// TO DO: exclude specific ICD10 codes from the analysis.
func initializeIcd10AnalysisMapsCCSR(icd10ToCssrMap map[string]ccsrCategory) (map[string][]int, map[int]Icd10Entry, int) {
	analysisIdMap := map[string][]int{}      // maps icd 10 code to analysis IDs
	analysisIcd10Map := map[int]Icd10Entry{} // maps analysis ID to a medical Name
	ccsrIDMap := map[string]int{}
	ctr := 0 //serves as analysis ID generator
	icd10ToExclude := getIcd10CodesToExcludeFromAnalysis()
	for icd10Code, ccsr := range icd10ToCssrMap {
		if _, ok := icd10ToExclude[icd10Code[0:1]]; ok {
			continue
		}
		var ids []int
		for id, name := range ccsr.categories {
			var ccsrID int
			var ok bool
			if ccsrID, ok = ccsrIDMap[id]; !ok {
				ccsrID = ctr
				analysisIcd10Map[ctr] = Icd10Entry{Name: name}
				ccsrIDMap[id] = ccsrID
				ctr++
			}
			ids = append(ids, ccsrID)
		}
		analysisIdMap[icd10Code] = ids
	}
	extra := getNonICD10CodesToAddToAnalysis()
	for code, name := range extra {
		analysisIcd10Map[ctr] = Icd10Entry{Name: name}
		analysisIdMap[code] = []int{ctr}
		ctr++
	}
	fmt.Println("Mapped ", len(icd10ToCssrMap), " ICD10 codes to ", ctr, " analysis IDs")
	return analysisIdMap, analysisIcd10Map, ctr
}

type icd10AnalysisMapsFromCCSR struct {
	Icd10Map          map[int]Icd10Entry // map analysis DID -> Icd10Entry
	NofDiagnosisCodes int                // nr of different diagnosis codes
	DIDMap            map[string][]int   // maps ICD10 Code onto multiple DIDs
}

type icd10AnalysisMapsFromXML struct {
	Icd10Map          map[int]Icd10Entry // map analysis DID -> Icd10Entry
	NofDiagnosisCodes int                // nr of different diagnosis codes
	DIDMap            map[string]int     // map ICD10 Code -> DID
}

func (analysisMap icd10AnalysisMapsFromXML) getDID(icd10Name string) int {
	if v, ok := analysisMap.DIDMap[icd10Name]; ok {
		return v
	} else {
		return -1
	}
}

func (analysisMap icd10AnalysisMapsFromCCSR) getDID(icd10DID string) []int {
	if v, ok := analysisMap.DIDMap[icd10DID]; ok {
		return v
	}
	return nil
}

func (analysisMap icd10AnalysisMapsFromXML) GetICDCode(did int) string {
	for icd10Code, didCode := range analysisMap.DIDMap {
		if didCode == did {
			return icd10Code
		}
	}
	return ""
}

func (analysisMap icd10AnalysisMapsFromCCSR) GetICDCode(did int) string {
	for icd10Code, didCodes := range analysisMap.DIDMap {
		for _, didCode := range didCodes {
			if didCode == did {
				return icd10Code
			}
		}
	}
	return ""
}

func (analysisMap icd10AnalysisMapsFromXML) getIdMap() map[int]string {
	res := map[int]string{}
	for icd10Code, didCode := range analysisMap.DIDMap {
		res[didCode] = icd10Code
	}
	return res
}

func (analysisMap icd10AnalysisMapsFromCCSR) getIdMap() map[int]string {
	res := map[int]string{}
	for icd10Code, didCodes := range analysisMap.DIDMap {
		for _, didCode := range didCodes {
			res[didCode] = icd10Code
		}
	}
	return res
}

// AnalysisMaps represent maps extracted from the input that map analysis IDs onto medical terms and vice versa. This is
// an interface that defines several methods. getICDCode returns for a did the original id in the input for the
// diagnostic event. fillInPatientDiagnoses creates for a given diagnosis identifier from the input a Diagnosis object
// and adds it to a patient's list of diagnoses.
type AnalysisMaps interface {
	fillInPatientDiagnoses(patient *Patient, DidString string, date DiagnosisDate) int
	fillInNonICDPatientDiagnoses(patient *Patient, infoMap map[string]*TreatmentInfo) int
	GetICDCode(did int) string
	getIdMap() map[int]string
}

func (analysisMap icd10AnalysisMapsFromXML) fillInPatientDiagnoses(patient *Patient, DIDString string, date DiagnosisDate) int {
	DID := analysisMap.getDID(DIDString)
	if DID == -1 {
		return 1 // icd10 diagnosis excluded from analysis
	}
	diagnosis := &Diagnosis{PID: patient.PID, DID: DID, Date: date}
	patient.AddDiagnosis(diagnosis)
	return 0
}

func (analysisMap icd10AnalysisMapsFromCCSR) fillInPatientDiagnoses(patient *Patient, DIDString string, date DiagnosisDate) int {
	DIDs := analysisMap.getDID(DIDString)
	if DIDs == nil {
		return 1 // icd10 code excluded from analysis
	}
	for _, DID := range DIDs {
		diagnosis := &Diagnosis{PID: patient.PID, DID: DID, Date: date}
		patient.AddDiagnosis(diagnosis)
	}
	return 0
}

func (analysisMap icd10AnalysisMapsFromXML) fillInNonICDPatientDiagnoses(patient *Patient, infoMap map[string]*TreatmentInfo) int {
	nonIcd := 0
	if info, ok := infoMap[patient.PIDString]; ok {
		if info.RCDate != nil {
			diagnosis := &Diagnosis{PID: patient.PID, DID: analysisMap.DIDMap["C98"], Date: *info.RCDate}
			nonIcd = 1
			patient.AddDiagnosis(diagnosis)
		}
		if info.MVACDate != nil {
			nonIcd = 1
			diagnosis := &Diagnosis{PID: patient.PID, DID: analysisMap.DIDMap["C99"], Date: *info.MVACDate}
			patient.AddDiagnosis(diagnosis)
		}
		if info.IVTDate != nil {
			nonIcd = 1
			diagnosis := &Diagnosis{PID: patient.PID, DID: analysisMap.DIDMap["C100"], Date: *info.IVTDate}
			patient.AddDiagnosis(diagnosis)
		}
	}
	return nonIcd
}

func (analysisMap icd10AnalysisMapsFromCCSR) fillInNonICDPatientDiagnoses(patient *Patient, infoMap map[string]*TreatmentInfo) int {
	nonIcd := 0
	if info, ok := infoMap[patient.PIDString]; ok {
		if info.RCDate != nil {
			dids := analysisMap.DIDMap["C98"]
			for _, did := range dids {
				diagnosis := &Diagnosis{PID: patient.PID, DID: did, Date: *info.RCDate}
				nonIcd = 1
				patient.AddDiagnosis(diagnosis)
			}
		}
		if info.MVACDate != nil {
			dids := analysisMap.DIDMap["C99"]
			for _, did := range dids {
				diagnosis := &Diagnosis{PID: patient.PID, DID: did, Date: *info.MVACDate}
				nonIcd = 1
				patient.AddDiagnosis(diagnosis)
			}
		}
		if info.IVTDate != nil {
			dids := analysisMap.DIDMap["C100"]
			for _, did := range dids {
				nonIcd = 1
				diagnosis := &Diagnosis{PID: patient.PID, DID: did, Date: *info.IVTDate}
				patient.AddDiagnosis(diagnosis)
			}
		}
	}
	return nonIcd
}

// initializeIcd10AnalysisMaps returns a map ICD10 DID -> internal analysis DID and a map analysis DID ->
// medical Name for an ICD10 Hierarchy passed as xml file and a requested hierarchy Level.
func initializeIcd10AnalysisMapsFromXML(file string, level int) icd10AnalysisMapsFromXML {
	icd10MapFromXml := initializeIcd10NameMap(file) // map ICD10 DID -> ICD 10 Name (medical desc, Categories, Level)
	analysisIdMap, icd10Map, ctr := initializeIcd10AnalysisMaps(icd10MapFromXml, level)
	return icd10AnalysisMapsFromXML{DIDMap: analysisIdMap, Icd10Map: icd10Map, NofDiagnosisCodes: ctr}
}

// initializeIcd10AnalysisMapsFromCCSR returns a map ICD10 -> []{internal analysis DID} and map analysis DID -> medical
// Name for ICD10 CCSR categorization passed as a csv file.
func initializeIcd10AnalysisMapsFromCCSR(file string) icd10AnalysisMapsFromCCSR {
	icd10ToCssrMap := initializeIcd10ToCCSRMap(file) // map ICD10 Code -> CCSR Name
	analysisIdMap, icd10Map, ctr := initializeIcd10AnalysisMapsCCSR(icd10ToCssrMap)
	return icd10AnalysisMapsFromCCSR{DIDMap: analysisIdMap, Icd10Map: icd10Map, NofDiagnosisCodes: ctr}
}

//Parsing patient information.

// parseTriNetXPatientData parses a file with patient information from the TriNetX database. Input: a patient file in csv
// format, a desired number of age groups to initialize cohorts. Diagnoses of the patient need to be filled in after
// parsing the diagnoses file.
func parseTriNetXPatientData(file string, nofCohortAges int) (*PatientMap, int) {
	//open file
	csvFile, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer func() {
		if err := csvFile.Close(); err != nil {
			panic(err)
		}
	}()
	patientMap := &PatientMap{PIDMap: map[int]*Patient{}, PIDStringMap: map[string]int{}}
	maxYOB := 1850
	minYOB := time.Now().Year() - 1
	deathCr := 0
	regions := map[string]int{} //counts per region
	regionIds := map[string]int{}
	//parse file
	reader := csv.NewReader(csvFile)
	//the header is omitted from the TriNetX file, but is should be: patient_id, sex, race, ethnicity, year_of_birth,
	//age_at_death, patient_regional_location, postal_code, marital_status, reason_yob_missing, month_year_death,
	//source_id
	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			panic(err)
		}
		var yob int
		if yob, err = strconv.Atoi(record[4]); err != nil {
			continue //skip patients without year of birth
		}
		pidString := record[0]
		patientMap.Ctr++      // avoid using 0 as PID
		pid := patientMap.Ctr //analysis ID
		var sex int
		if record[1] == "M" {
			sex = Male
			patientMap.MaleCtr++
		}
		if record[1] == "F" {
			sex = Female
			patientMap.FemaleCtr++
		}
		dateOfDeathString := record[10]
		var dateOfDeath *DiagnosisDate
		if len(dateOfDeathString) == 6 {
			year, err := strconv.Atoi(dateOfDeathString[0:4])
			if err == nil {
				month, err := strconv.Atoi(dateOfDeathString[4:6])
				if err == nil {
					deathCr++
					dateOfDeath = &DiagnosisDate{
						Year:  year,
						Month: month,
						Day:   1, //unknown, default to 1
					}
				}
			}
		}
		region := record[6]
		if _, ok := regions[region]; !ok {
			regions[region] = 0
			regionIds[region] = len(regionIds)
		} else {
			regions[region]++
		}
		patient := Patient{
			PID:       pid,
			PIDString: pidString,
			YOB:       yob,
			CohortAge: 0,
			Sex:       sex,
			Diagnoses: []*Diagnosis{},
			DeathDate: dateOfDeath,
			Region:    regionIds[region],
		}
		patientMap.PIDMap[pid] = &patient
		patientMap.PIDStringMap[pidString] = pid
		maxYOB = utils.MaxInt(yob, maxYOB)
		minYOB = utils.MinInt(yob, minYOB)
	}
	// initialize patient age groups
	ageRange := float64(maxYOB-minYOB) / float64(nofCohortAges)
	ageRange = math.Ceil(ageRange)
	if nofCohortAges > 1 {
		for _, p := range patientMap.PIDMap {
			p.CohortAge = int(math.Floor(float64(p.YOB-minYOB) / float64(ageRange)))
		}
	}
	fmt.Println("Parsed patient data.")
	fmt.Print("Parsed ", patientMap.Ctr, " patients with year of birth known ")
	fmt.Print("of which ", patientMap.FemaleCtr, " females and ")
	fmt.Println(patientMap.MaleCtr, "males; and of which ", deathCr, " have a known date of death.")
	fmt.Println("Year of birth oldest patient:", minYOB)
	fmt.Println("Year of birth youngest patient:", maxYOB)
	fmt.Println("Patients are of ", len(regions), " regions: ")
	for region, nr := range regions {
		fmt.Print(region, ": ", nr, ", ")
	}
	fmt.Println("")
	return patientMap, len(regions)
}

//Parsing patient diagnoses

// parseTriNetXDiagnosisDate turns a TriNetX date string into DiagnosisDate object.
func parseTriNetXDiagnosisDate(date string) DiagnosisDate {
	year, err := strconv.Atoi(date[0:4])
	if err != nil {
		panic(err)
	}
	month, err := strconv.Atoi(date[5:7])
	if err != nil {
		panic(err)
	}
	day, err := strconv.Atoi(date[8:10])
	if err != nil {
		panic(err)
	}
	return DiagnosisDate{Year: year, Month: month, Day: day}
}

// TriNetXEventOfInterest checks if the ICD10 code is related to bladder cancer
func TriNetXEventOfInterest(icd10ID string) bool {
	if icd10ID == "Z85.1" {
		return true
	}
	if len(icd10ID) >= 3 && icd10ID[0:3] == "C67" {
		return true
	}
	return false
}

// TreatmentInfo implements a structure for storing the dates of certain bladder cancer treatments.
type TreatmentInfo struct {
	RCDate   *DiagnosisDate //Date of radical cystectomy
	MVACDate *DiagnosisDate //Date of MVAC chemotherapy
	IVTDate  *DiagnosisDate //Date of intravesical therapy
}

// parseTriNetXTreatmentFile parses a csv file that contains information of patient's treatments at different time stamps.
// It returns a map from PID -> TreatmentInfo.
func parseTriNetXTreatmentFile(fileName string) map[string]*TreatmentInfo {
	result := map[string]*TreatmentInfo{}
	file, err := os.Open(fileName)
	if err != nil {
		panic(err)
	}
	defer func() {
		if err := file.Close(); err != nil {
			panic(err)
		}
	}()
	reader := csv.NewReader(file)
	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			panic(err)
		}
		PIDString := record[0]
		var rcDate, mvacDate, ivtDate *DiagnosisDate
		if len(record[10]) == 10 { // valid date
			d := parseTriNetXDiagnosisDate(record[10])
			rcDate = &d
		}
		if len(record[11]) == 10 {
			d := parseTriNetXDiagnosisDate(record[11])
			mvacDate = &d
		}
		if len(record[13]) == 10 {
			d := parseTriNetXDiagnosisDate(record[13])
			rcDate = &d
		}
		result[PIDString] = &TreatmentInfo{RCDate: rcDate, MVACDate: mvacDate, IVTDate: ivtDate}
	}
	return result
}

// parseTrinetXPatientDiagnoses parses a csv file containing patient diagnoses. It fills in those diagnoses for the given
// patients. It uses the icd10AnalysisMap to assign internal analysis DID to the diagnoses.
// TO DO: Handle ICD09 diagnoses.
func parseTrinetXPatientDiagnoses(diagnosesFile, treatmentInfoFile string, patients *PatientMap, icd10AnalysisMap AnalysisMaps, icd9ToIcd10Map map[string]string) {
	file, err := os.Open(diagnosesFile)
	if err != nil {
		panic(err)
	}
	defer func() {
		if err := file.Close(); err != nil {
			panic(err)
		}
	}()
	reader := csv.NewReader(file)
	ctr := 0 //for counting the number of parsed diagnoses
	ctrID09 := 0
	ctrExcl := 0
	EOICtr := 0
	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			panic(err)
		}
		ctr++
		PIDString := record[0]
		patient, ok := GetPatient(PIDString, patients)
		if !ok {
			continue //skip unknown patients
		}
		DIDCodeSystem := record[2]
		DIDString := record[3]
		if DIDCodeSystem != "ICD-10-CM" {
			// try to remap ICD9 code to ICD10 codes
			if DIDString, ok = icd9ToIcd10Map[DIDString]; !ok {
				continue // skip unkown ICD9 codes
			}
			ctrID09++
		}
		date := parseTriNetXDiagnosisDate(record[7])

		nr := icd10AnalysisMap.fillInPatientDiagnoses(patient, DIDString, date)
		if nr > 0 {
			ctrExcl++
			continue
		}
		//Check if diagnosis is event of interest.
		if patient.EOIDate == nil && TriNetXEventOfInterest(DIDString) {
			EOICtr++
			patient.EOIDate = &date // mark first event of interest (e.g. bladder cancers diagnosis)
		}
	}
	var nonICD10DiagnosesMap map[string]*TreatmentInfo
	nonICDCtr := 0
	if treatmentInfoFile != "" {
		nonICD10DiagnosesMap = parseTriNetXTreatmentFile(treatmentInfoFile)
		for _, patient := range patients.PIDMap {
			//fill in non ICD10 diagnoses derived from procedure info
			r := icd10AnalysisMap.fillInNonICDPatientDiagnoses(patient, nonICD10DiagnosesMap)
			nonICDCtr = nonICDCtr + r
		}
	}
	for _, patient := range patients.PIDMap {
		SortDiagnoses(patient)
		CompactDiagnoses(patient)
	}
	fmt.Println("Parsed diagnosis data.")
	fmt.Print("Parsed ", ctr, " diagnoses ")
	fmt.Println("of which ", ctrID09, " ICD09 diagnoses and ", ctr-ctrID09, " ICD10 diagnoses, and ", ctrExcl, " diagnoses excluded from analysis")
	fmt.Println("and of which ", EOICtr, " events of interest.")
	fmt.Println("Parsed non ICD diagnoses for: ", nonICDCtr, " patients.")
}

func ParseTriNetXData(name, patientFile, diagnosisFile, diagnosisInfoFile, treatmentInfoFile string, nofCohortAges,
	level int, minYears, maxYears float64, icd9ToIcd10File string, filters []PatientFilter) (*Experiment, *PatientMap) {
	// parse data
	// fill in patients
	patients, nofRegions := parseTriNetXPatientData(patientFile, nofCohortAges)
	// fill in icd10 to analysis map
	var analysisMaps AnalysisMaps
	var nofDiagnosisCodes int
	var icd10Map map[int]Icd10Entry
	var idMap map[int]string
	if filepath.Ext(diagnosisInfoFile) == ".xml" {
		maps := initializeIcd10AnalysisMapsFromXML(diagnosisInfoFile, level)
		analysisMaps = maps
		nofDiagnosisCodes = maps.NofDiagnosisCodes
		icd10Map = maps.Icd10Map
		idMap = maps.getIdMap()
	}
	if filepath.Ext(diagnosisInfoFile) == ".csv" || filepath.Ext(diagnosisInfoFile) == ".CSV" {
		maps := initializeIcd10AnalysisMapsFromCCSR(diagnosisInfoFile)
		analysisMaps = maps
		nofDiagnosisCodes = maps.NofDiagnosisCodes
		icd10Map = maps.Icd10Map
		idMap = maps.getIdMap()
	}
	icd9ToIcd10Map := map[string]string{}
	if icd9ToIcd10File != "" {
		icd9ToIcd10Map = parseIcd9ToIcd10Mapping(icd9ToIcd10File)
	}
	// fill in diagnoses for patients
	parseTrinetXPatientDiagnoses(diagnosisFile, treatmentInfoFile, patients, analysisMaps, icd9ToIcd10Map)
	// Apply patient filter
	patients = ApplyPatientFilters(filters, patients)
	fmt.Println("Filtered down to: ", len(patients.PIDMap), " patients.")
	// create cohorts
	cohorts := InitCohorts(patients, nofCohortAges, nofRegions, nofDiagnosisCodes)
	mergedCohort := MergeCohorts(cohorts)
	exp := Experiment{
		NofAgeGroups:      nofCohortAges,
		Level:             level,
		NofDiagnosisCodes: nofDiagnosisCodes,
		DxDRR:             MakeDxDRR(nofDiagnosisCodes),
		DxDPatients:       MakeDxDPatients(nofDiagnosisCodes),
		DPatients:         mergedCohort.DPatients,
		Cohorts:           cohorts,
		Name:              name,
		Icd10Map:          icd10Map,
		NofRegions:        nofRegions,
		IdMap:             idMap,
		FCtr:              patients.FemaleCtr,
		MCtr:              patients.MaleCtr,
	}
	return &exp, patients
}

// opening json file with ICD09 -> ICD10 mapping

func parseIcd9ToIcd10Mapping(file string) map[string]string {
	jsonFile, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer jsonFile.Close()
	fmt.Println("Parsing ICD9 to ICD10 mapping from a json file.")
	jsonBytes, _ := io.ReadAll(jsonFile)
	var mapping map[string]string
	json.Unmarshal(jsonBytes, &mapping)
	return mapping
}

// TumorInfo is a struct for storing bladder cancer tumor information concerning: tumor size, tumor lymph nodes, tumor
// metastasis
type TumorInfo struct {
	TStage, NStage, MStage, Stage string
	Date                          DiagnosisDate
}

// getTumorStage converts tumor size, number of lymph nodes, and metastatis Level into an overall cancer stage.
// T stages: Ta,T1,Tis,T2,T3,T4
// N stages: N0,N1,N2,N3
// M stages: M0,M1
// Stage 0a: Ta,N0,M0
// Stage 0is:Tis,N0,M0 known as carcinoma in situ (CIS)
// Stage I: T1,N0,M0
// Stage II: T2,N0,M0
// Stage IIIA: T3a,T3b, or T4a,N0,M0 --or-- T1 to T4a,N1,M0
// Stage IIIB: T1 to T4a, N2 or N3, M0
// Stage IVA: T4b,any N,M0 or any T, any N, M1a
// Stage IVB: any T, any N, M1b
func getTumorStage(tStage, nStage, mStage string) string {
	if nStage == "N0" && mStage == "M0" {
		switch tStage {
		case "Ta":
			return "0a"
		case "Tis":
			return "0is"
		case "T1":
			return "I"
		case "T2":
			return "II"
		case "T3a", "T3b", "T4a":
			return "IIIA"
		}
	}
	if nStage == "N1" && mStage == "M0" {
		switch tStage {
		case "T1", "T1a", "T1b", "T2", "T2a", "T2b", "T3", "T3a", "T3b", "T4a":
			return "IIIA"
		}
	}
	if (nStage == "N2" || nStage == "N3") && mStage == "M0" {
		switch tStage {
		case "T1", "T1a", "T1b", "T2", "T2a", "T2b", "T3", "T3a", "T3b", "T4", "T4a":
			return "IIIB"
		}
	}
	if tStage == "T4b" && mStage == "M0" {
		return "IVA"
	}
	if mStage == "M1a" {
		return "IVA"
	}
	if mStage == "M1b" {
		return "IVB"
	}
	return tStage + nStage + mStage
}

// tumorIsCISStage checks if tumor is flat or carcinoma in situ (CIS).
func tumorIsCISStage(tumor *TumorInfo) bool {
	return tumor.Stage == "0is"
}

// ParsetTriNetXTumorData parses the tumor data from a csv file and returns a map PIDString -> []*TumorInfo.
func ParsetTriNetXTumorData(fileName string) map[string][]*TumorInfo {
	file, err := os.Open(fileName)
	if err != nil {
		panic(err)
	}
	defer func() {
		if err := file.Close(); err != nil {
			panic(err)
		}
	}()
	result := map[string][]*TumorInfo{}
	reader := csv.NewReader(file)
	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			panic(err)
		}
		tumorSite := strings.Split(record[4], ".")
		if tumorSite[0] == "C67" { //only record bladder cancer information
			PIDString := record[0]
			date := parseTriNetXDiagnosisDate(record[1])
			tumorSizeInfo := strings.Split(record[10], "_")
			numberOfLymphNodesInfo := strings.Split(record[11], "_")
			metastaticInfo := strings.Split(record[12], "_")
			if len(tumorSizeInfo) == 1 || len(numberOfLymphNodesInfo) == 1 || len(metastaticInfo) == 1 {
				continue
			}
			tumor := &TumorInfo{Date: date, TStage: tumorSizeInfo[1], NStage: numberOfLymphNodesInfo[1],
				MStage: metastaticInfo[1]}
			tumor.Stage = getTumorStage(tumorSizeInfo[1], numberOfLymphNodesInfo[1], metastaticInfo[1])
			if ts, ok := result[PIDString]; ok {
				result[PIDString] = append(ts, tumor)
			} else {
				result[PIDString] = []*TumorInfo{tumor}
			}
		}
	}
	printTumorInfoSummary(result)
	return result
}

func printTumorInfoSummary(tumorInfo map[string][]*TumorInfo) {
	fmt.Println("Parsed tumor info. Found tumor info for: ", len(tumorInfo), " patients.")
	ctr := map[string]int{}
	for _, tumors := range tumorInfo {
		for _, tumor := range tumors {
			ctr[tumor.TStage]++
			ctr[tumor.NStage]++
			ctr[tumor.MStage]++
		}
	}
	stages := []string{}
	for stage, _ := range ctr {
		stages = append(stages, stage)
	}
	sort.Strings(stages)
	for _, stage := range stages {
		fmt.Println("For stage: ", stage, ": ", ctr[stage], " entries.")
	}
}

func printTumorInfo(tumorInfo map[int][]*TumorInfo) {
	fmt.Println("Nr of patients with tumor info: ", len(tumorInfo))
	for pid, infos := range tumorInfo {
		for _, info := range infos {
			fmt.Println("Patient with pid: ", pid, " has TStage: ", info.TStage, " NStage: ", info.NStage, " MStage: ",
				info.MStage, " Global stage: ", info.Stage)
		}
	}
}
