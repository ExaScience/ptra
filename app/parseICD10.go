package app

import (
	"encoding/csv"
	"encoding/xml"
	"fmt"
	"io"
	"io/ioutil"
	"os"
	"ptra/trajectory"
	"slices"
)

//Package ptra implements a patient trajectory analysis tool.
//The ptra program has 3 data inputs:
//A file mapping diagnosis ID (DID) -> medical name.
//A file with patient info, associating a patient ID (PID) with YOB, sex, etc.
//A file with diagnosis info, mapping PID -> DID, date.

//TriNetX Database stores diagnoses data using a mix of ICD10 and ICD9 codes.
//We have an additional file that maps ICD9 IDs -> ICD10 IDs.
//We can download the ICD10 ID -> medical name from https://www.cms.gov/medicare/icd-10/2022-icd-10-cm as an xml file.
//TriNetX stores patient info as a csv file, as well as the diagnoses info.
//For the medical name mapping, we can also use the ICD10 -> CCSR Mapping which maps ICD10 onto 530 categories with medical meaning. This
//mapping can be downloaded from https://www.hcup-us.ahrq.gov/toolssoftware/ccsr/dxccsr.jsp#download as a CSV file. This
//mapping performs a mapping ICD10 -> CCSR categpries -> medical name.

//Parsing ICD10 names hierarchy from xml
//Structs for unmarshalling ICD10 xml data
//Structure of an ICD10 code: ABC.XYZ(D): up to 7 characters
//The xml file structures these codes using:
//<chapter> <section> <desc> </desc> </section> </chapter>
//with chapter the first level of the diagnosis code (A), section the second level (D) and desc the rest (C.XYZ(D))

//The ptra program needs a names map that maps DID -> medical name. The following code extracts a name map from an ICD10
//hierarchy and a given level.

// CDC ICD10 hierarchy
// parseIcd10HierarchyFromXML parses the xml file with the ICD10 hierarchy into an icd10Hierarchy object.
func parseIcd10HierarchyFromXml(file string) icd10Hierarchy {
	fmt.Println("Parsing ICD10 code hierarchy from XML file: ", file)
	//open file
	xmlFile, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer xmlFile.Close()
	xmlFileBytes, _ := ioutil.ReadAll(xmlFile)
	//unmarshall
	icd10Hierarchy := icd10Hierarchy{}
	xml.Unmarshal(xmlFileBytes, &icd10Hierarchy)
	return icd10Hierarchy
}

// printIcd10Hierarchy prints an ICD10 hierarchy parsed from an XML file.
func printIcd10Hierarchy(hierarchy icd10Hierarchy) {
	fmt.Println("Printing ICD10 code hierarchy.")
	// count # DID per level
	ctr1, ctr2, ctr3, ctr4, ctr5, ctr6, ctr7 := 0, 0, 0, 0, 0, 0, 0
	for _, chap := range hierarchy.Chapters {
		// level 1
		ctr1++
		fmt.Println("Chapter: ", chap.Desc)
		for _, section := range chap.Sections {
			// level 2
			ctr2++
			fmt.Println("Section: ", section.Desc)
			for _, diag := range section.Diagnoses {
				// level 3
				ctr3++
				fmt.Println(diag.Name, " : ", diag.Desc)
				if len(diag.Diagnoses) == 0 {
					continue
				}
				for _, diag := range diag.Diagnoses {
					// level 4
					ctr4++
					fmt.Println(diag.Name, " : ", diag.Desc)
					if len(diag.Diagnoses) == 0 {
						continue
					}
					for _, diag := range diag.Diagnoses {
						// level 5
						ctr5++
						fmt.Println(diag.Name, " : ", diag.Desc)
						if len(diag.Diagnoses) == 0 {
							continue
						}
						for _, diag := range diag.Diagnoses {
							// level 6
							ctr6++
							fmt.Println(diag.Name, " : ", diag.Desc)
							if len(diag.Diagnoses) == 0 {
								continue
							}
							// level 7
							ctr7++
							fmt.Println(diag.Name, " : ", diag.Desc)
						}
					}
				}
			}
		}
	}
	fmt.Println("#ICD10 codes/descriptors per level: ")
	fmt.Println("Lvl 0: ", ctr1, " Lvl 1: ", ctr2, " Lvl 2: ",
		ctr3, " Lvl 3: ", ctr4, " Lvl 4: ", ctr5, " Lvl 5: ", ctr6, " Lvl 6: ", ctr7)
}

type cdcICD10HierarchyXML struct {
	XMLName xml.Name `xml:"ICD10CM.tabular"`
}

// diag captures the lowest levels of the ICD10 code
type diag struct {
	Name      string `xml:"name"` //Unique diagnosis ID (DID) in ICD10 encoding
	Desc      string `xml:"desc"` //A medical name/description for a DID
	Diagnoses []diag `xml:"diag"` //A diagnosis can be split into more detailed diagnoses descriptors.
}

// section captures the second level of the ICD10 code
type section struct {
	Desc      string `xml:"desc"`    //A medical name/description for a DID.
	Id        string `xml:"id,attr"` //Unique diagnosis ID (DID) in ICD10 encoding
	Diagnoses []diag `xml:"diag"`    //A diagnosis can be split into more detailed diagnoses descriptors.
}

// chapter captures the first (highest) level of the ICD10 code
type chapter struct {
	XmlName  xml.Name  `xml:"chapter"`
	Desc     string    `xml:"desc"`
	Sections []section `xml:"section"`
}

// icd10Hierarchy contains the full xml table with the ICD10 code hierarchy.
type icd10Hierarchy struct {
	XmlName  xml.Name  `xml:"ICD10CM.tabular"`
	Chapters []chapter `xml:"chapter"`
}

// WHO ICD10
type SuperClass struct {
	XMLName xml.Name `xml:"SuperClass"`
	Code    string   `xml:"code,attr"`
}

type whoRubric struct {
	XMLName xml.Name `xml:"Rubric"`
	Label   string   `xml:"Label"`
}

type whoClass struct {
	XmlName    xml.Name    `xml:"Class"`
	Kind       string      `xml:"kind,attr"`
	Code       string      `xml:"code,attr"`
	Rubrics    []whoRubric `xml:"Rubric"`
	SuperClass SuperClass  `xml:"SuperClass"`
}

func (class whoClass) isChapter() bool {
	return class.Kind == "chapter"
}

func (class whoClass) isBlock() bool {
	return class.Kind == "block"
}

func (class whoClass) isCategory() bool {
	return class.Kind == "category"
}

type whoIcd10Hierarchy struct {
	XmlName  xml.Name   `xml:"ClaML"` //not being parsed??
	Chapters []whoClass `xml:"Class"`
}

type whoICD10HierarchyXML struct {
	XmlName xml.Name `xml:"Class"`
}

func ParseWhoIcd10HierarchyFromXml(file string) whoIcd10Hierarchy {
	fmt.Println("Parsing who ICD10 code hierarchy from XML file: ", file)
	//open file
	xmlFile, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer xmlFile.Close()
	xmlFileBytes, _ := ioutil.ReadAll(xmlFile)
	//unmarshall
	whoICD10Hierarchy := whoIcd10Hierarchy{}
	xml.Unmarshal(xmlFileBytes, &whoICD10Hierarchy)
	return whoICD10Hierarchy
}

// CheckIcd10HierarchyXMLFile checks if an xml file contains the CDC ICD10 hierarchy or the WHO ICD10 hierarchy
// Returns a string "cdc" or "who" respectively
func CheckIcd10HierarchyXMLFile(file string) string {
	//open file
	xmlFile, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer xmlFile.Close()
	xmlFileBytes, _ := io.ReadAll(xmlFile)
	cdcHierarchy := cdcICD10HierarchyXML{}
	xml.Unmarshal(xmlFileBytes, &cdcHierarchy)
	if cdcHierarchy.XMLName.Local == "ICD10CM.tabular" {
		return "cdc"
	}
	whoHierarchy := whoICD10HierarchyXML{}
	xml.Unmarshal(xmlFileBytes, &whoHierarchy)
	if whoHierarchy.XmlName.Local == "Class" {
		return "who"
	}
	fmt.Println(whoHierarchy)
	return "unknown"
}

// ICD10 structures
// icd10Name is a struct for containing a medical name + level + the categories of a DID in ICD10 encoding.
type icd10Name struct {
	name       string    //medical name for a DID in ICD10 encoding
	categories [6]string //the names of the ICD10 encoding higher and lower in the hierarchy.
	level      int       //the ICD10 hierarchy level of this name.
}

type icd10Table struct {
	namesMap map[string]icd10Name //maps ICD10 DID to a medical name, level, and categories to which it belongs.
}

// selectParentCategory returns the ICD10 name of the parent category to which a DID belongs.
func selectParentCategory(name icd10Name) string {
	if name.level > 0 {
		return name.categories[name.level-1]
	}
	return "None"
}

func printIcd10NameMap(table map[string]icd10Name) {
	fmt.Println("ICD10 Name map: ")
	for id, name := range table {
		fmt.Println(id, " : ", name)
	}
}

// initializeIcd10NameMap initializes a name map for ICD10 DID -> medical name, level, and categories it belongs to.
func initializeIcd10NameMap(file string) map[string]icd10Name {
	target := CheckIcd10HierarchyXMLFile(file)
	if target == "who" {
		return InitializeWhoIcd10NameMap(file)
	}
	if target == "cdc" {
		return initializeCdcIcd10NameMap(file)
	}
	panic("Unknown file format ICD10 hierarchy")
}

func initializeCdcIcd10NameMap(file string) map[string]icd10Name {
	icd10NameMap := map[string]icd10Name{} //maps ICD10 DID to a medical name, level, and categories to which it belongs.
	icd10Hierarchy := parseIcd10HierarchyFromXml(file)
	for _, chap := range icd10Hierarchy.Chapters {
		category0 := chap.Desc
		for _, section := range chap.Sections {
			category1 := section.Desc
			// manually unrolled loop since we know hierarchy is max 7 levels, otherwise recursive code
			for _, diag := range section.Diagnoses {
				if len(diag.Diagnoses) == 0 {
					icd10Name := icd10Name{name: diag.Desc,
						categories: [6]string{category0, category1, "NONE", "NONE", "NONE", "NONE"}, level: 2}
					icd10NameMap[diag.Name] = icd10Name
					continue
				}
				category2 := diag.Desc
				for _, diag := range diag.Diagnoses {
					if len(diag.Diagnoses) == 0 {
						icd10Name := icd10Name{name: diag.Desc,
							categories: [6]string{category0, category1, category2, "NONE", "NONE", "NONE"},
							level:      3}
						icd10NameMap[diag.Name] = icd10Name
						continue
					}
					category3 := diag.Desc
					for _, diag := range diag.Diagnoses {
						if len(diag.Diagnoses) == 0 {
							ICD10Name := icd10Name{name: diag.Desc,
								categories: [6]string{category0, category1, category2, category3, "NONE", "NONE"},
								level:      4}
							icd10NameMap[diag.Name] = ICD10Name
							continue
						}
						category4 := diag.Desc
						for _, diag := range diag.Diagnoses {
							if len(diag.Diagnoses) == 0 {
								ICD10Name := icd10Name{name: diag.Desc,
									categories: [6]string{category0, category1, category2, category3, category4, "NONE"},
									level:      5}
								icd10NameMap[diag.Name] = ICD10Name
								continue
							}
							category5 := diag.Desc
							for _, diag := range diag.Diagnoses {
								ICD10Name := icd10Name{name: diag.Desc,
									categories: [6]string{category0, category1, category2, category3, category4, category5},
									level:      6}
								icd10NameMap[diag.Name] = ICD10Name
							}
						}
					}
				}
			}
		}
	}
	return icd10NameMap
}

func InitializeWhoIcd10NameMap(file string) map[string]icd10Name {
	icd10NameMap := map[string]icd10Name{} //maps ICD10 DID to a medical name, level, and categories to which it belongs.
	whoIcd10Hierarchy := ParseWhoIcd10HierarchyFromXml(file)
	//map code to chapter
	chapters := map[string]whoClass{}
	blocks := map[string]whoClass{}
	categories := map[string]whoClass{}

	for _, class := range whoIcd10Hierarchy.Chapters {
		if class.isChapter() {
			chapters[class.Code] = class
			continue
		}
		if class.isBlock() {
			blocks[class.Code] = class
			continue
		}
		if class.isCategory() {
			categories[class.Code] = class
		}
	}

	getParentCategory := func(id string) (whoClass, bool) {
		v, ok := categories[id]
		if ok {
			return v, true
		}
		v, ok = blocks[id]
		if ok {
			return v, true
		}
		v, ok = chapters[id]
		if ok {
			return v, true
		}
		return whoClass{}, false
	}

	for _, cat := range categories {
		parent, ok := getParentCategory(cat.SuperClass.Code)
		pcats := []whoClass{parent}
		//grab category parents
		for ok {
			parent, ok = getParentCategory(parent.SuperClass.Code)
			if ok {
				pcats = append(pcats, parent)
			}
		}
		entry := [6]string{"NONE", "NONE", "NONE", "NONE", "NONE", "NONE"}
		pctr := 0
		slices.Reverse(pcats)
		for _, pcat := range pcats {
			entry[pctr] = pcat.Rubrics[0].Label
			pctr++
		}
		icd10Name := icd10Name{name: cat.Rubrics[0].Label, categories: entry, level: pctr}
		icd10NameMap[cat.Code] = icd10Name
	}
	return icd10NameMap
}

// initializeIcd10AnalysisIDMap creates a map ICD10 DID -> analysis DID and a map analysis ID -> medical name. This is
// useful to remap diagnosis codes used in the input to a higher level in the ICD10 hierarchy. E.g "typhoid fever" and
// "cholera" are both "infectuous intestinal diseases", so they could both be identified as such during the analysis.
// This can be interesting to obtain more global patient trajectories/clusters.
func intializeIcd10AnalysisMaps(icd10NameMap map[string]icd10Name, level int,
	icd10ToExclude map[string]bool) (map[string]int, map[int]string, int) {
	analysisIdMap := map[string]int{}       // maps icd 10 code to analysis ID
	analysisNameMap := map[int]string{}     // maps analysis ID to a medical name
	nameToAnalysisIdMap := map[string]int{} // maps medical name to analysis ID
	ctr := 0                                //serves as analysis ID generator
	for icd10Code, icd10Name := range icd10NameMap {
		if _, ok := icd10ToExclude[icd10Name.categories[0]]; ok {
			// code to exclude from analysis
			continue
		}
		var name string
		if level == icd10Name.level || level > icd10Name.level {
			name = icd10Name.name
		} else {
			name = icd10Name.categories[level]
		}
		// may already have seen name, because of level
		newID, ok := nameToAnalysisIdMap[name]
		if !ok {
			newID = ctr
			ctr++
			analysisNameMap[newID] = name
			nameToAnalysisIdMap[name] = newID
		}
		analysisIdMap[icd10Code] = newID
	}
	extra := getNonICD10CodesToAddToAnalysis()
	for code, name := range extra {
		analysisNameMap[ctr] = name
		nameToAnalysisIdMap[name] = ctr
		analysisIdMap[code] = ctr
		ctr++
	}
	fmt.Println("Mapped ", len(icd10NameMap), " ICD10 codes to ", ctr, " analysis IDs of level ", level)
	return analysisIdMap, analysisNameMap, ctr
}

// ccsrCategory is a struct for containing CCSR categories, encoding medically meaningful names for a DID in ICD10
// encoding.
type ccsrCategory struct {
	name       string            //default CCSR category/medical name
	id         string            //CCSR ID for default category
	categories map[string]string //Up to 6 different CCSR categories an ICD10 code is mapped to
}

type icd10ToCCSRTable map[string]ccsrCategory //maps ICD10 DID to its CCSR categories

// ccsrIcd10ToProperIcd10 transforms the ICD10 code from a ccsr file into a proper ICD10 code. The ICD10 codes in the
// ccsr file are stored without the ".", so this needs to be added to be able to compare to any other data that uses
// ICD10 codes. Also removes superfluous quotes stored in the ccsr file.
func ccsrIcd10ToProperIcd10(code string) string {
	return code[1:4] + "." + code[4:len(code)-1]
}

// initializeIcd10NameMapFromCCSR initializes a name map for ICD10 DID -> CCSR categories (medical names)
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
		//fill in unique CSSR alternative categories, up to 6 possible
		for i := 6; i <= 17; i = i + 2 {
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

// initializeIcd10AnalysisMapsCCSR creates a map ICD10 DID -> [analysis DID] and a map analysis ID -> medical name,
// starting from a CCSR mapping, which maps ICD10 codes onto medical meaningful categories.
// Each icd10 code can be mapped to multiple ccsr categories, and therefore to multiple analysis IDs.
// TO DO: exclude specific ICD10 codes from the analysis.
func initializeIcd10AnalysisMapsCCSR(icd10ToCssrMap map[string]ccsrCategory,
	icd10ToExclude map[string]bool) (map[string][]int, map[int]string, int) {
	analysisIdMap := map[string][]int{} // maps icd 10 code to analysis IDs
	analysisNameMap := map[int]string{} // maps analysis ID to a medical name
	ccsrIDMap := map[string]int{}
	ctr := 0 //serves as analysis ID generator
	for icd10Code, ccsr := range icd10ToCssrMap {
		if _, ok := icd10ToExclude[icd10Code[0:1]]; ok {
			continue
		}
		ids := []int{}
		for id, name := range ccsr.categories {
			var ccsrID int
			var ok bool
			if ccsrID, ok = ccsrIDMap[id]; !ok {
				ccsrID = ctr
				analysisNameMap[ctr] = name
				ccsrIDMap[id] = ccsrID
				ctr++
			}
			ids = append(ids, ccsrID)
		}
		analysisIdMap[icd10Code] = ids
	}
	extra := getNonICD10CodesToAddToAnalysis()
	for code, name := range extra {
		analysisNameMap[ctr] = name
		analysisIdMap[code] = []int{ctr}
		ctr++
	}
	fmt.Println("Mapped ", len(icd10ToCssrMap), " ICD10 codes to ", ctr, " analysis IDs")
	return analysisIdMap, analysisNameMap, ctr
}

type icd10AnalysisMapsFromCCSR struct {
	NameMap           map[int]string   // map analysis DID -> medical name
	NofDiagnosisCodes int              // nr of different diagnosis codes
	DIDMap            map[string][]int // maps ICD10 Code onto multiple DIDs
}

type icd10AnalysisMapsFromXML struct {
	NameMap           map[int]string // map analysis DID -> medical name
	NofDiagnosisCodes int            // nr of different diagnosis codes
	DIDMap            map[string]int // map ICD10 Code -> DID
}

func (analysisMap icd10AnalysisMapsFromXML) getDID(icd10DID string) int {
	if v, ok := analysisMap.DIDMap[icd10DID]; ok {
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
	fillInPatientDiagnoses(patient *trajectory.Patient, DidString string, date trajectory.DiagnosisDate) int
	fillInNonICDPatientDiagnoses(patient *trajectory.Patient, infoMap map[string]*TreatmentInfo) int
	GetICDCode(did int) string
	getIdMap() map[int]string
}

func (analysisMap icd10AnalysisMapsFromXML) fillInPatientDiagnoses(patient *trajectory.Patient, DIDString string,
	date trajectory.DiagnosisDate) int {
	DID := analysisMap.getDID(DIDString)
	if DID == -1 {
		return 1 // icd10 diagnosis excluded from analysis
	}
	diagnosis := &trajectory.Diagnosis{PID: patient.PID, DID: DID, Date: date}
	trajectory.AddDiagnosis(patient, diagnosis)
	return 0
}

func (analysisMap icd10AnalysisMapsFromCCSR) fillInPatientDiagnoses(patient *trajectory.Patient, DIDString string,
	date trajectory.DiagnosisDate) int {
	DIDs := analysisMap.getDID(DIDString)
	if DIDs == nil {
		return 1 // icd10 code excluded from analysis
	}
	for _, DID := range DIDs {
		diagnosis := &trajectory.Diagnosis{PID: patient.PID, DID: DID, Date: date}
		trajectory.AddDiagnosis(patient, diagnosis)
	}
	return 0
}

func (analysisMap icd10AnalysisMapsFromXML) fillInNonICDPatientDiagnoses(patient *trajectory.Patient,
	infoMap map[string]*TreatmentInfo) int {
	nonIcd := 0
	if info, ok := infoMap[patient.PIDString]; ok {
		if info.RCDate != nil {
			diagnosis := &trajectory.Diagnosis{PID: patient.PID, DID: analysisMap.DIDMap["C98"], Date: *info.RCDate}
			nonIcd = 1
			trajectory.AddDiagnosis(patient, diagnosis)
		}
		if info.MVACDate != nil {
			nonIcd = 1
			diagnosis := &trajectory.Diagnosis{PID: patient.PID, DID: analysisMap.DIDMap["C99"], Date: *info.MVACDate}
			trajectory.AddDiagnosis(patient, diagnosis)
		}
		if info.IVTDate != nil {
			nonIcd = 1
			diagnosis := &trajectory.Diagnosis{PID: patient.PID, DID: analysisMap.DIDMap["C100"], Date: *info.IVTDate}
			trajectory.AddDiagnosis(patient, diagnosis)
		}
	}
	return nonIcd
}

func (analysisMap icd10AnalysisMapsFromCCSR) fillInNonICDPatientDiagnoses(patient *trajectory.Patient,
	infoMap map[string]*TreatmentInfo) int {
	nonIcd := 0
	if info, ok := infoMap[patient.PIDString]; ok {
		if info.RCDate != nil {
			dids := analysisMap.DIDMap["C98"]
			for _, did := range dids {
				diagnosis := &trajectory.Diagnosis{PID: patient.PID, DID: did, Date: *info.RCDate}
				nonIcd = 1
				trajectory.AddDiagnosis(patient, diagnosis)
			}
		}
		if info.MVACDate != nil {
			dids := analysisMap.DIDMap["C99"]
			for _, did := range dids {
				diagnosis := &trajectory.Diagnosis{PID: patient.PID, DID: did, Date: *info.MVACDate}
				nonIcd = 1
				trajectory.AddDiagnosis(patient, diagnosis)
			}
		}
		if info.IVTDate != nil {
			dids := analysisMap.DIDMap["C100"]
			for _, did := range dids {
				nonIcd = 1
				diagnosis := &trajectory.Diagnosis{PID: patient.PID, DID: did, Date: *info.IVTDate}
				trajectory.AddDiagnosis(patient, diagnosis)
			}
		}
	}
	return nonIcd
}

// initializeIcd10AnalysisMaps returns a map ICD10 DID -> internal analysis DID and a map analysis DID ->
// medical name for an ICD10 Hierarchy passed as xml file and a requested hierarchy level.
func initializeIcd10AnalysisMapsFromXML(file string, level int, icd10ToExclude map[string]bool) icd10AnalysisMapsFromXML {
	icd10NameMapFromXml := initializeIcd10NameMap(file) // map ICD10 DID -> ICD 10 Name (medical desc, categories, level)
	analysisIdMap, analysisNameMap, ctr := intializeIcd10AnalysisMaps(icd10NameMapFromXml, level, icd10ToExclude)
	return icd10AnalysisMapsFromXML{DIDMap: analysisIdMap, NameMap: analysisNameMap, NofDiagnosisCodes: ctr}
}

// initializeIcd10AnalysisMapsFromCCSR returns a map ICD10 -> []{internal analysis DID} and map analysis DID -> medical
// name for ICD10 CCSR categorization passed as a csv file.
func initializeIcd10AnalysisMapsFromCCSR(file string, icd10ToExclude map[string]bool) icd10AnalysisMapsFromCCSR {
	icd10ToCssrMap := initializeIcd10ToCCSRMap(file) // map ICD10 Code -> CCSR Name
	analysisIdMap, analysisNameMap, ctr := initializeIcd10AnalysisMapsCCSR(icd10ToCssrMap, icd10ToExclude)
	return icd10AnalysisMapsFromCCSR{DIDMap: analysisIdMap, NameMap: analysisNameMap, NofDiagnosisCodes: ctr}
}
