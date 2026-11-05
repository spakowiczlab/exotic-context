#!/usr/bin/env python3
"""Merge species-level Bracken reports into counts and non-host fractions.

The final fraction table uses the HMP vaginal reference layout. Each sample's
denominator is the total Bracken-estimated reads over retained species. This is
not an absolute microbial abundance estimate, and unclassified reads are not in
the denominator. The input reports and excluded taxa are retained separately.
"""

import argparse
import csv
import math
import os
import re
import sys
from array import array
from pathlib import Path


def load_samples(path):
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not {"sample", "status"} <= set(reader.fieldnames):
            raise ValueError("sample status file must contain sample and status columns")
        rows = list(reader)
    samples = [row["sample"] for row in rows if row["status"] == "done"]
    ids = [row["sample"] for row in rows]
    if not samples or len(set(ids)) != len(ids):
        raise ValueError("missing or duplicate samples in status file")
    if any(not re.fullmatch(r"SRS[0-9]+", sid) for sid in ids):
        raise ValueError("unexpected sample identifier")
    if any(row["status"] not in {"done", "excluded"} for row in rows):
        raise ValueError("unfinished or failed samples in status file")
    return sorted(samples)


def load_report(path, sample_idx, n_samples, taxa):
    excluded = {"9606": 0, "32630": 0}
    kept = 0
    seen = set()
    required = {"name", "taxonomy_id", "taxonomy_lvl", "new_est_reads"}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not required <= set(reader.fieldnames):
            raise ValueError(f"unexpected Bracken report columns: {path}")
        for row in reader:
            tid, name, level = row["taxonomy_id"].strip(), row["name"].strip(), row["taxonomy_lvl"].strip()
            if not re.fullmatch(r"[0-9]+", tid) or not name or level != "S":
                raise ValueError(f"invalid species row in {path}: {tid} {name} {level}")
            if tid in seen:
                raise ValueError(f"duplicate taxid {tid} in {path}")
            seen.add(tid)
            try:
                count = int(row["new_est_reads"])
            except (TypeError, ValueError) as exc:
                raise ValueError(f"invalid estimated count for taxid {tid} in {path}") from exc
            if count < 0:
                raise ValueError(f"negative estimated count for taxid {tid} in {path}")
            if tid in excluded:
                excluded[tid] += count
                continue
            if tid not in taxa:
                taxa[tid] = (name, level, array("q", [0]) * n_samples)
            old_name, old_level, counts = taxa[tid]
            if name != old_name or level != old_level:
                raise ValueError(f"inconsistent taxonomy for {tid}: {old_name!r} vs {name!r}")
            counts[sample_idx] = count
            kept += count
    if kept <= 0:
        raise ValueError(f"no retained estimated species reads in {path}")
    return kept, excluded


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-status", type=Path, required=True)
    parser.add_argument("--reports-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--suffix", default=".bracken.tsv", help="report filename suffix")
    parser.add_argument("--expected-samples", type=int, default=466)
    args = parser.parse_args()
    if "/" in args.suffix or not args.suffix:
        raise ValueError("suffix must be a filename suffix, without slashes")
    samples = load_samples(args.sample_status)
    if len(samples) != args.expected_samples:
        raise ValueError(f"expected {args.expected_samples} included samples, found {len(samples)}")
    reports = [args.reports_dir / f"{sid}{args.suffix}" for sid in samples]
    missing = [str(path) for path in reports if not path.is_file() or path.stat().st_size == 0]
    if missing:
        raise ValueError(f"missing/empty Bracken reports: {len(missing)}; examples: {missing[:5]}")
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    names = ["HMP_gut_species_counts_nonhost.tsv", "HMP_gut_species_relative_abundance_nonhost_renormalized.tsv", "HMP_gut_species_sample_QC.tsv"]
    if any((out / name).exists() for name in names):
        raise ValueError("output file already exists; select a new output directory")
    taxa = {}
    totals, exclusions = [], []
    for idx, report in enumerate(reports):
        kept, excluded = load_report(report, idx, len(samples), taxa)
        totals.append(kept)
        exclusions.append(excluded)
    if not taxa or any(not math.isfinite(x) or x <= 0 for x in totals):
        raise ValueError("no valid taxa or sample denominator")

    paths = [out / name for name in names]
    temps = [out / f".{name}.{os.getpid()}.partial" for name in names]
    for p in temps:
        if p.exists():
            raise ValueError(f"temporary file already exists: {p}")
    try:
        with temps[0].open("x", newline="", encoding="utf-8") as cf, temps[1].open("x", newline="", encoding="utf-8") as rf:
            count_out = csv.writer(cf, delimiter="\t", lineterminator="\n")
            fraction_out = csv.writer(rf, delimiter="\t", lineterminator="\n")
            header = ["taxonomy_id", "name", "taxonomy_lvl", *samples]
            count_out.writerow(header)
            fraction_out.writerow(header)
            for tid in sorted(taxa, key=int):
                name, level, counts = taxa[tid]
                count_out.writerow([tid, name, level, *counts])
                fraction_out.writerow([tid, name, level, *(format(c / total, ".12g") for c, total in zip(counts, totals))])
        with temps[2].open("x", newline="", encoding="utf-8") as qf:
            qc = csv.writer(qf, delimiter="\t", lineterminator="\n")
            qc.writerow(["sample", "retained_estimated_reads", "human_9606_estimated_reads", "synthetic_32630_estimated_reads"])
            for sample, total, ex in zip(samples, totals, exclusions):
                qc.writerow([sample, total, ex["9606"], ex["32630"]])
        for temp, dest in zip(temps, paths):
            if dest.exists():
                raise ValueError(f"output appeared during writing: {dest}")
            os.rename(temp, dest)
    except BaseException:
        for temp in temps:
            temp.unlink(missing_ok=True)
        raise
    print(f"Samples: {len(samples)}; retained species: {len(taxa)}")
    for path in paths:
        print(path)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
