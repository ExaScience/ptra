# ptra\_mimic\_iv — MIMIC-IV Use Case

Standalone `ptra` application for trajectory analysis on [MIMIC-IV](https://physionet.org/content/mimiciv/) hospital data. Two event sources are supported, selected via `--source`:

- **`hcpcs`** (default) — HCPCS procedure events from `hcpcsevents.csv`
- **`diagnoses_icd`** — ICD diagnosis events from `diagnoses_icd.csv`, with ICD-9 → ICD-10 conversion

For general background on the algorithm and the analysis pipeline, see the [main README](../README.md).

## Input Data

`patients.csv` is always required. The remaining files depend on `--source`:

| `--source` | Required files | Purpose |
|------------|---------------|---------|
| `hcpcs` | `hcpcsevents.csv` | HCPCS procedure events with `chartdate` as the event date |
| `diagnoses_icd` | `diagnoses_icd.csv`, `admissions.csv`, `d_icd_diagnoses.csv` | ICD diagnoses joined with admission dates; names from the definitions table |

| File | Used columns |
|------|-------------|
| `patients.csv` | `subject_id`, `gender`, `anchor_age`, `anchor_year` |
| `hcpcsevents.csv` | `subject_id`, `chartdate`, `hcpcs_cd`, `short_description` |
| `diagnoses_icd.csv` | `subject_id`, `hadm_id`, `icd_code`, `icd_version` |
| `admissions.csv` | `hadm_id`, `admittime` |
| `d_icd_diagnoses.csv` | `icd_code`, `long_title` |

Columns are resolved by header name, so column order does not matter.

### ICD-9 → ICD-10 conversion

When using `--source diagnoses_icd`, rows with `icd_version=9` are converted to ICD-10 using a JSON mapping file passed via `--ICD9ToICD10File` (same format as the TriNetX `--ICD9ToICD10File`). Rows that cannot be converted are skipped. If no mapping file is provided, all ICD-9 rows are skipped.

## Build

```bash
go build -o bin/ptra_mimic_iv ./ptra_mimic_iv/
```

## Usage

```
ptra_mimic_iv <hospDir> <outputPath> [flags]
```

### Examples

```bash
# HCPCS events (default)
./bin/ptra_mimic_iv /path/to/mimiciv/hosp/ ./output/

# ICD diagnoses with ICD-9 → ICD-10 conversion
./bin/ptra_mimic_iv /path/to/mimiciv/hosp/ ./output/ \
  --source diagnoses_icd --ICD9ToICD10File ICD_9_to_10.json
```

### Flags

| Flag | Default | Description |
|------|---------|-------------|
| `--source` | hcpcs | Event source: `hcpcs` or `diagnoses_icd` |
| `--ICD9ToICD10File` | | JSON mapping ICD-9 → ICD-10 (for `diagnoses_icd` source) |
| `--nofAgeGroups` | 5 | Number of age groups for cohort stratification |
| `--minYears` | 0.5 | Minimum years between consecutive events in a pair |
| `--maxYears` | 5.0 | Maximum years between consecutive events in a pair |
| `--minPatients` | 5 | Minimum patients for last event in a trajectory |
| `--minTrajectoryLength` | 3 | Minimum number of events in a trajectory |
| `--maxTrajectoryLength` | 5 | Maximum number of events in a trajectory |
| `--iter` | 400 | Monte Carlo iterations for RR calculation |
| `--RR` | 1.0 | Minimum relative-risk score for pairs |
| `--saveRR` | | Save RR matrix to file (reusable with `--loadRR`) |
| `--loadRR` | | Load RR matrix from a previous run |
| `--cluster` | false | Run MCL clustering on resulting trajectories |
| `--mclPath` | | Path to the `mcl` binary (required when `--cluster` is set) |
| `--pfilters` | id | Patient filters: `male`, `female`, `age70+`, `age70-` (comma-separated) |
| `--nrOfThreads` | 0 | Number of threads (0 = Go default) |

## Differences from the TriNetX Use Case

- **Two event sources** — HCPCS codes or ICD diagnoses, selected via `--source`. TriNetX always uses ICD codes through an XML hierarchy with level selection; here each code is its own event type.
- **ICD-9 → ICD-10 conversion** — uses the same JSON format as TriNetX (`--ICD9ToICD10File`). The lookup handles both dotted and undotted code formats since MIMIC-IV stores codes without dots.
- **Diagnosis dates from admissions** — `diagnoses_icd.csv` has no date column; the admission date (`admittime`) from `admissions.csv` is used via `hadm_id` join.
- **No tumor/treatment filters** — the TriNetX-specific cancer-stage filters are not available; only demographic filters apply.
- **Demographics from anchor years** — MIMIC-IV de-identifies dates; year of birth is estimated from `anchor_year − anchor_age`.
