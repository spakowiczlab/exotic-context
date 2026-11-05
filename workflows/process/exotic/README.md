# exotic: HMP gut processing

HMP gut processing for the [exotic pipeline](https://github.com/spakowiczlab/exotic), using KrakenUniq/Bracken. See the upstream repository for pipeline documentation, package usage, and database information. The exact upstream revision used for these results remains to be recorded in [`profiles/exotic-001.yaml`](../../../profiles/exotic-001.yaml).

Analysis code originally uploaded under `HMP_gut_analysis_code/`, relocated here without changing the submitted scripts or supporting metadata.

| Path | Purpose |
| --- | --- |
| `hmp_gut_pilot_submit.sh` | Submit a pilot KrakenUniq/Bracken analysis |
| `hmp_gut_kraken_full_submit.py` | Validate prior inputs/pilot state and submit full processing and merge jobs |
| `scripts/worker.sbatch` | Submitted Slurm worker script |
| `scripts/merge.sbatch` | Submitted Slurm merge job |
| `scripts/merge.py` | Merge per-sample outputs |
| `manifest.tsv` | Original sample/read-count manifest snapshot |
| `tables/sample_status.tsv` | Original prior-processing status snapshot |
| `task_ids.txt` | Original array task sample IDs |

The metadata snapshots remain with the code to preserve the submitted bundle and its internal layout. They are not the standardized `datasets/hmp1/samples.tsv` or `runs.tsv` tables, and their status values refer to the original processing context.

These scripts retain their original HPC paths, Slurm settings, sample counts, and pilot job assumptions. Relocation does not make them portable or verify their provenance. The submitters read existing cluster metadata and create runtime scripts and files under an external run directory; moving this repository directory does not change those paths. Inspect and configure the scripts before execution. No jobs were submitted as part of this reorganization.
