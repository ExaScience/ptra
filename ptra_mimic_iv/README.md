# ptra\_mimic\_iv — MIMIC-IV Use Case

Standalone `ptra` application for trajectory analysis on [MIMIC-IV](https://physionet.org/content/mimiciv/) hospital data. It uses **HCPCS procedure events** (`hcpcsevents.csv`) as the temporal events for building patient trajectories.

For general background on the algorithm and the analysis pipeline, see the [main README](../README.md).

## Input Data

Two CSV files from the MIMIC-IV `hosp` module are required:

| File              | Used columns | Purpose |
|-------------------|-------------|---------|
| `patients.csv`    | `subject_id`, `gender`, `anchor_age`, `anchor_year` | Patient demographics; YOB is estimated as `anchor_year − anchor_age` |
| `hcpcsevents.csv` | `subject_id`, `chartdate`, `hcpcs_cd`, `short_description` | Temporal events; each unique HCPCS code becomes an event type |

Columns are resolved by header name, so column order does not matter.

## Build

```bash
go build -o bin/ptra_mimic_iv ./ptra_mimic_iv/
```

## Usage

```
ptra_mimic_iv <hospDir> <outputPath> [flags]
```

`hospDir` is the path to the MIMIC-IV `hosp` directory that contains `patients.csv` and `hcpcsevents.csv`. The program locates both files automatically.

### Example

```bash
./bin/ptra_mimic_iv \
  /path/to/mimiciv/hosp/ \
  ./output/ \
  --nofAgeGroups 10 --minPatients 50 --iter 400 --saveRR output/rr.csv
```

### Flags

| Flag | Default | Description |
|------|---------|-------------|
| `--nofAgeGroups` | 6 | Number of age groups for cohort stratification |
| `--minYears` | 0.5 | Minimum years between consecutive events in a pair |
| `--maxYears` | 5.0 | Maximum years between consecutive events in a pair |
| `--minPatients` | 100 | Minimum patients for last event in a trajectory |
| `--minTrajectoryLength` | 3 | Minimum number of events in a trajectory |
| `--maxTrajectoryLength` | 5 | Maximum number of events in a trajectory |
| `--iter` | 10000 | Monte Carlo iterations for RR calculation |
| `--RR` | 1.0 | Minimum relative-risk score for pairs |
| `--saveRR` | | Save RR matrix to file (reusable with `--loadRR`) |
| `--loadRR` | | Load RR matrix from a previous run |
| `--cluster` | false | Run MCL clustering on resulting trajectories |
| `--mclPath` | | Path to the `mcl` binary (required when `--cluster` is set) |
| `--pfilters` | id | Patient filters: `male`, `female`, `age70+`, `age70-` (comma-separated) |
| `--nrOfThreads` | 0 | Number of threads (0 = Go default) |

## Differences from the TriNetX Use Case

- **No ICD code hierarchy or level selection** — HCPCS codes are used directly; each unique code is its own event type with `short_description` as the human-readable name.
- **No ICD-9 → ICD-10 conversion** — not applicable to HCPCS.
- **No tumor/treatment filters** — the TriNetX-specific cancer-stage filters (`NMIBC`, `MIBC`, etc.) are not available; only demographic filters apply.
- **Demographics from anchor years** — MIMIC-IV de-identifies dates; year of birth is estimated from `anchor_year − anchor_age`.

