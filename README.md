# micro-context

Public microbiome datasets reprocessed with versioned lab workflows and reference databases for consistent comparison.

This repository catalogs source data, curated metadata, processing profiles, and release provenance. Pipelines remain in their own repositories; wrappers here will orchestrate their use. Consistent processing reduces computational differences, but does not remove differences introduced by sample collection, extraction, library preparation, sequencing, or study design. EXOTIC and bioBakery outputs are separate processing tracks.

## Start here

- [Prepare or add a dataset](docs/adding-datasets.md)
- [Define a processing profile](docs/processing-profiles.md)
- [Find and use outputs](docs/using-the-data.md)
- [HMP starter dataset](datasets/hmp1/README.md)

## Organization

```text
datasets/                 Source descriptions, samples, sequencing runs, metadata mappings
profiles/                 Pinned pipeline, software, reference database, and parameter records
workflows/{download,process,export}/  Future orchestration wrappers
schemas/                  JSON Schema contracts for metadata and release records
releases/<dataset>/<profile>/        Versioned release manifests linking external data
scripts/                  Scaffold validation
tests/fixtures/          Small synthetic metadata examples
```

## Current status

Existing HMP count/abundance tables are currently at the repository root, with submitted analysis code under `HMP_gut_analysis_code/`. These files are preserved as uploaded and have not yet been audited or linked to a verified processing release.

The HMP records under `datasets/` are draft templates for cataloging existing EXOTIC and bioBakery processing. They contain no real samples or verified processing versions. `hmp1` is a provisional dataset ID: confirm the source phase and subset before populating it. No download or processing commands are implemented yet.

Large data live outside Git. Store FASTQs, reference databases, intermediates, and large result tables in durable storage, and record URLs, sizes, and checksums in manifests. Small curated metadata belong here. Local staging directories `data/`, `results/`, and `scratch/` are ignored.

## Validation

Run `python3 scripts/validate_metadata.py` (standard library only). This checks table headers, unique identifiers, sample/run relationships, and the synthetic fixture. JSON Schemas describe row objects after parsing TSV files; YAML manifests can be checked against the corresponding schemas with a YAML-aware JSON Schema validator.

## Releases and attribution

A published release identifies the dataset selection, processing profile, exact source commit, and external output files. Never overwrite a published release; create a new release manifest and retain the old one. Draft manifests are incomplete and must not be treated as downloadable releases.

Code is MIT licensed; source datasets retain their original terms. Cite the original datasets, pipelines, databases, and the specific micro-context release used. See [CITATION.cff](CITATION.cff).
