# Processing profiles

A profile defines the computational conditions under which outputs can be compared. Pin pipeline revision, individual tool versions, container digests or environment locks, all reference databases and checksums, parameters, preprocessing, multi-run handling, and export conventions. Store exact invocation/configuration and run logs with release provenance.

Profile IDs (`exotic-001`, `biobakery-001`) are stable identifiers, not software versions. Drafts may be edited. Once used in a published release, changing code, databases, parameters, or output conventions requires a new profile ID. Apply the same profile across datasets intended for comparison.

Keep EXOTIC and bioBakery tracks separate. Cross-pipeline feature mappings are explicit versioned derivatives, with unmapped features reported. Do not assume that matching taxonomic names imply identical feature definitions or that abundance units are interchangeable.

For bioBakery, resolve the reported “v4” into the exact workflow, constituent tools, and databases used. Do not infer component versions from that label.

Pipeline implementation stays upstream. Add thin wrappers under `workflows/` only after the actual pipeline entry points and execution environment are known.
