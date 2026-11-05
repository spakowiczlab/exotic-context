# HMP starter dataset (draft)

This directory is reserved for the HMP subset already processed by the lab through exotic and bioBakery. Confirm that `hmp1` is the correct phase before adding records; rename the dataset directory and references if necessary.

## Prepare the existing data

1. Complete `dataset.yaml`: source accessions, citations, selection criteria, source terms, and metadata retrieval date.
2. Add one row per biological sample to `samples.tsv`. Preserve source IDs and subject IDs; use `unknown` for unavailable values. Scope subject IDs to this dataset.
3. Add one row per sequencing run to `runs.tsv`. Encode multiple runs of one sample as separate rows. Paired reads occupy the two URL/checksum columns of the same row. Document any different file structure before extending the table.
4. Record original-to-curated metadata mappings in `metadata_mapping.tsv`; preserve original metadata in external storage with provenance.
5. Recover exact processing details in `profiles/exotic-001.yaml` and `profiles/biobakery-001.yaml`. The bioBakery version hint is user-reported and unverified.
6. Inventory existing outputs in the draft release manifests. Include intended samples that failed or were excluded, with reasons.

The TSVs contain headers only. No HMP accessions, sample counts, database releases, or QC thresholds have been assumed.

Source discovery: https://www.hmpdacc.org/ and https://www.ncbi.nlm.nih.gov/sra
