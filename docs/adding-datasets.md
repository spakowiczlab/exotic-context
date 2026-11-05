# Adding datasets and preparing existing results

Use a stable dataset ID that names the study and relevant subset. Copy the HMP starter directory, update IDs, and write explicit inclusion/exclusion criteria. Record source accessions and citations, retrieval date, source terms, and original metadata locations.

## Metadata conventions

`samples.tsv` has one biological sample per row. `sample_id` is unique within the dataset; qualify it with `dataset_id` when combining studies. Preserve original sample IDs. `subject_id` is dataset-scoped and links repeated samples; never infer identity across studies. `collection_time` uses ISO 8601 when an actual timestamp is available; describe relative time and units in the metadata mapping instead of inventing dates. Missing metadata use `unknown`.

`runs.tsv` has one sequencing run per row, linked to `sample_id`. `layout` is `single` or `paired`. Use retrievable URLs and SHA-256 checksums for the actual input files. For single-end runs, set both read-2 fields to `unknown`. For interleaved reads or other file structures, extend and document the schema first. Record original archive checksums and transformations in supplemental provenance when they differ from the input files.

Map source fields to standardized fields in `metadata_mapping.tsv`, recording transformations, units, and source missing-value conventions. Extra columns are allowed, but document their meaning. Preserve source metadata separately.

## Existing processed datasets

Inventory the files before rerunning anything. Recover pipeline versions, database releases/checksums, settings, logs, sample selection, and merge behavior. Mark unrecoverable provenance `unknown`; do not substitute current defaults. A legacy product with incomplete provenance may be useful but must retain that limitation in its release notes.

Stage large files in ignored `data/` or `results/`, then upload to durable storage. Add URLs, sizes, and SHA-256 values to release outputs. Review `git status` before committing; ignored file patterns do not replace a file inventory.

Run `python3 scripts/validate_metadata.py`. Before publishing, verify source access, sample counts, all intended sample outcomes, output checksums, profile completeness, feature conventions, and dataset attribution. A release cannot be complete while required provenance is unknown.
