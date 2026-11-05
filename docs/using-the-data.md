# Using releases

Choose a dataset and processing profile, then find a published release manifest under `releases/<dataset>/<profile>/`. Draft manifests are placeholders. A release should link sample metadata, feature tables, feature annotations, per-sample QC, processing status, and provenance. Each external file needs a retrievable URL, byte size, SHA-256 checksum, format, and units where relevant.

Before combining datasets, check profile identity, assay, body site, study design, repeated subjects, QC, feature identifiers, table orientation, and units. Computational harmonization does not remove study-level technical or biological differences. Filter by documented criteria and retain excluded/failed sample records.

Preserve native pipeline outputs. Standardized exports must document orientation, taxonomic rank or functional namespace, units, filtering, zero versus missing values, and any feature mapping. Record those transformations in provenance.

## Publishing

Retain draft status until inputs, profile, provenance, and outputs are verified. Set `status: published`, record the exact source commit and UTC creation timestamp, and attach a release tag. For later changes add `manifest-r002.yaml` (and subsequent versions) with a new release ID; retain previous manifests and external files unchanged. New processing conditions require a new profile directory.

Cite the original study/data source, pipelines, databases, and specific micro-context release. The repository MIT license does not change source-data terms.
