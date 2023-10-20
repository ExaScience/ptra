#!/bin/sh

# Look for the patient file
cd input || exit

# Determine patient file
if [ -e "$PATIENT_FILE" ]; then
    echo "Patient file $PATIENT_FILE found."
elif [ -n "$PATIENT_FILE" ]; then
   echo "Patient file '$PATIENT_FILE' not found."
   exit 1
elif [ -e "patients.csv" ]; then
   PATIENT_FILE="patients.csv"
else
  echo "patients.csv NOT FOUND"
  exit 1
fi

# Determine diagnostics file
if [ -e "$DIAGNOSIS_INFO_FILE" ]; then
    echo "Diagnosis file $DIAGNOSIS_INFO_FILE found."
elif [ -n "$DIAGNOSIS_INFO_FILE" ]; then
   echo "Diagnosis file $DIAGNOSIS_INFO_FILE not found."
   exit 1
elif [ -e "diagnosisInfo.xml" ]; then
   DIAGNOSIS_INFO_FILE="diagnosisInfo.xml"
elif [ -e "diagnosisInfo.csv" ]; then
  DIAGNOSIS_INFO_FILE="diagnosisInfo.csv"
else
  echo "diagnosisInfo.xml or diagnosisInfo.csv NOT FOUND"
  exit 1
fi

# Determine diagnoses file
if [ -e "$DIAGNOSES_FILE" ]; then
    echo "Diagnoses file $DIAGNOSES_FILE found."
elif [ -n "$DIAGNOSES_FILE" ]; then
   echo "Diagnoses file '$DIAGNOSES_FILE' not found."
   exit 1
elif [ -e "diagnoses.csv" ]; then
   DIAGNOSES_FILE="diagnoses.csv"
else
  echo "diagnoses.csv NOT FOUND"
  exit 1
fi

# Add flags
FLAGS=""
addFlag() {
  if [ -n "$1" ]; then
    FLAGS="$FLAGS --$2 $1"
  fi
}

addFlag "$NUMBER_OF_AGE_GROUPS" "nofAgeGroups"
addFlag "$LEVEL" "lvl"
addFlag "$MIN_PATIENTS" "minPatients"
addFlag "$MAX_YEARS" "maxYears"
addFlag "$MIN_YEARS" "minYears"
addFlag "$MAX_TRAJECTORY_LENGTH" "maxTrajectoryLength"
addFlag "$MIN_TRAJECTORY_LENGTH" "minTrajectoryLength"
addFlag "$NAME" "name"
addFlag "$ICD9_TO_ICD10_FILE" "ICD9ToICD10File"
addFlag "$CLUSTER" "cluster"
addFlag "$MCL_PATH" "mclPath"
addFlag "$ITER" "iter"
addFlag "$SAVE_RR" "saveRR"
addFlag "$LOAD_RR" "loadRR"
addFlag "$PFILTERS" "pfilters"
addFlag "$TUMOR_INFO" "tumorInfo"
addFlag "$TFILTERS" "tfilters"
addFlag "$TREATMENT_INFO" "treatmentInfo"
addFlag "$CLUSTER_GRANULARITIES" "clusterGranularities"
addFlag "$NUMBER_OF_THREADS" "nrOfThreads"
addFlag "$RR" "RR"

# Trim the flags
FLAGS=$(echo "$FLAGS" | sed 's/ *$//g')
FLAGS=$(echo "$FLAGS" | sed 's/--cluster 1/--cluster/g') # "--cluster" is a flag without parameter: to enable it, set it to "1"
echo "*$FLAGS*"
cd ..

# Input and output folder
if [ -n "$OUTPUT_FOLDER" ]; then
    echo "Using output: $OUTPUT_FOLDER"
    output_dir="$OUTPUT_FOLDER"
else
    output_dir="/output/"
fi
if [ -n "$INPUT_FOLDER" ]; then
    input_dir="$INPUT_FOLDER"
else
    input_dir="/input/"
fi

iso8601_datetime=$(date -u +"%Y%m%dT%H%M%SZ")
output_dir="$output_dir/$iso8601_datetime"

echo "Using input: $input_dir"
echo "Using output: $output_dir"

./ptra "$input_dir""$PATIENT_FILE" "$input_dir""$DIAGNOSIS_INFO_FILE" "$input_dir""$DIAGNOSES_FILE" "$output_dir" $FLAGS
