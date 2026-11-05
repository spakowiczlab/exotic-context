#!/usr/bin/env python3
"""Validate draft metadata tables; not a release readiness or URL check."""
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def read_table(path, kind):
    schema = json.loads((ROOT / "schemas" / f"{kind}.schema.json").read_text())
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = set(schema["required"]) - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path}: missing columns {sorted(missing)}")
        rows = list(reader)
    for line, row in enumerate(rows, 2):
        if None in row or any(not row.get(key) for key in schema["required"]):
            raise ValueError(f"{path}:{line}: malformed row or empty required value")
    return rows

def validate(directory):
    samples = read_table(directory / "samples.tsv", "samples")
    runs = read_table(directory / "runs.tsv", "runs")
    ids = [row["sample_id"] for row in samples]
    run_ids = [row["run_id"] for row in runs]
    if len(ids) != len(set(ids)) or len(run_ids) != len(set(run_ids)):
        raise ValueError(f"{directory}: duplicate sample or run ID")
    for row in runs:
        if row["sample_id"] not in ids:
            raise ValueError(f"{directory}: run references an absent sample")
        if row["layout"] not in {"single", "paired"}:
            raise ValueError(f"{directory}: layout must be single or paired")
        for key in ["read_1_sha256", "read_2_sha256"]:
            value = row[key]
            if value != "unknown" and (len(value) != 64 or any(c not in "0123456789abcdefABCDEF" for c in value)):
                raise ValueError(f"{directory}: invalid SHA-256")
        if row["layout"] == "single" and any(row[k] != "unknown" for k in ["read_2_url", "read_2_sha256"]):
            raise ValueError(f"{directory}: single-end run has read-2 values")
        if row["layout"] == "paired" and row["read_2_url"] == "unknown":
            raise ValueError(f"{directory}: paired run lacks a read-2 URL")
    print(f"OK {directory.relative_to(ROOT)}: {len(samples)} samples, {len(runs)} runs")

if __name__ == "__main__":
    for directory in sorted((ROOT / "datasets").iterdir()):
        if directory.is_dir():
            validate(directory)
    validate(ROOT / "tests" / "fixtures")
