# micro-context <img src="man/figures/micro-context-hex.png" align="right" width="150" alt="micro-context hex sticker" />

**Public microbiome data, processed consistently for the lab.**

micro-context is a shared collection of public microbiome datasets reprocessed with the Spakowicz Lab’s pipelines and reference databases. It gives lab members reusable reference data for putting new results in context, without repeating the work of finding inputs, reconstructing metadata, and processing the same studies for each project.

The project brings together dataset descriptions, analysis code, processing configurations, and links to reusable outputs. It is organized as a data resource and workflow repository; an installable R package is not currently provided.

## Goals

- **Make public data comparable.** Apply consistent software, reference databases, and parameters within each processing track to reduce computational differences between studies.
- **Process once and reuse.** Make curated metadata, abundance tables, and quality-control information available for recurring lab analyses.
- **Keep results traceable.** Record source samples and sequencing runs, pipeline revisions, database releases, parameters, and output checksums.
- **Make updates explicit.** Preserve released products and create new processing profiles when software, databases, or analytical conventions change.

Consistent processing does not eliminate differences caused by sample collection, extraction, library preparation, sequencing, or study design. Those factors remain relevant when selecting comparison datasets and interpreting results.

## Processing tracks and current data

The initial focus is Human Microbiome Project (HMP) data processed through [{exotic}](https://github.com/spakowiczlab/exotic) and [bioBakery workflows](https://github.com/biobakery/biobakery_workflows). Each track retains its own feature definitions and output conventions. Comparisons across studies should use the same processing profile; cross-pipeline mappings require their own documented transformations.

The submitted HMP {exotic} analysis code is in [`workflows/process/exotic/`](workflows/process/exotic/). It includes the original Slurm submission scripts and accompanying metadata snapshots. The scripts retain their original HPC settings and run-specific assumptions; see the workflow README before running them.

Existing HMP count and abundance tables remain at the repository root. They have not yet been linked to a verified release. The standardized metadata under `datasets/hmp1/`, the processing profiles, and the release manifests are draft templates awaiting the original run details. The reported bioBakery “v4” label still needs to be resolved into exact tool and database versions. Confirm the HMP phase and subset before finalizing the provisional `hmp1` dataset ID.

## Repository structure

```text
micro-context/
├── datasets/
│   └── hmp1/                   # Source description and curated metadata templates
│       ├── dataset.yaml        # Study scope, accessions, and selection criteria
│       ├── samples.tsv        # Biological samples and subject metadata
│       ├── runs.tsv           # Sequencing runs, input URLs, and checksums
│       └── metadata_mapping.tsv
├── profiles/                   # Pipeline, software, database, and parameter records
│   ├── exotic-001.yaml
│   └── biobakery-001.yaml
├── workflows/
│   ├── download/               # Reserved for input retrieval and verification
│   ├── process/
│   │   └── exotic/             # Submitted HMP {exotic} code and run snapshots
│   └── export/                 # Reserved for standardized output generation
├── releases/<dataset>/<profile>/ # Release manifests linking data and provenance
├── schemas/                    # JSON Schema contracts for metadata and releases
├── scripts/                    # Metadata validation utilities
├── tests/fixtures/             # Small synthetic metadata examples
├── docs/                       # Guides for preparing, processing, and using data
└── man/figures/                # README artwork
```

A **dataset** defines the source study and selected samples. A **processing profile** defines how those inputs are processed, including software and reference databases. A **release** ties the inputs and profile to specific output files and processing outcomes.

Pipeline development belongs in the upstream pipeline repositories. This repository holds dataset-specific execution code and provenance needed to reproduce and reuse the resulting data.

## Preparing and using data

Start with the [HMP preparation checklist](datasets/hmp1/README.md), or follow the guide to [add a dataset](docs/adding-datasets.md). Record exact processing details in the corresponding [processing profile](docs/processing-profiles.md), then inventory the outputs in a release manifest. Keep unknown provenance explicitly marked until it can be recovered.

For analysis, select a published release and follow the [data-use guide](docs/using-the-data.md). Release products should include sample metadata, feature tables and annotations, per-sample QC, processing outcomes, and provenance. Draft manifests are placeholders rather than published data releases.

Store FASTQs, reference databases, intermediates, and large result tables in durable external storage. Record their URLs, byte sizes, and checksums in manifests. Keep code, small curated metadata, and configurations in Git. Local staging directories `data/`, `results/`, and `scratch/` are ignored.

## Validation

```sh
python3 scripts/validate_metadata.py
```

The validator uses the Python standard library and checks required table headers, unique identifiers, sample/run relationships, and basic run-field consistency. It also checks a synthetic fixture. It does not verify external files, processing environments, or release readiness. The JSON Schemas describe parsed TSV rows and release manifests for additional validation.

## Releases and attribution

Published releases identify the dataset selection, processing profile, source commit, and external output files. Preserve published manifests and data; subsequent corrections or additions receive a new release ID. Changes to processing conditions require a new profile.

Code is [MIT licensed](LICENSE); source datasets retain their original terms. Cite the original datasets, pipelines, reference databases, and the specific micro-context release used. See [CITATION.cff](CITATION.cff).
