#!/usr/bin/env python3
"""Submit one gut KrakenUniq/Bracken array and a completeness-checked merge.

The array waits for the previously submitted pilot to succeed. By default it
allows up to 200 concurrent sample jobs; OSC may queue tasks under its limits.
Each sample retains KrakenUniq and Bracken reports. The final merged files use
estimated species reads and 0-1 nonhost relative abundances.
"""

import argparse
import csv
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile


MERGE_SOURCE = r'''#!/usr/bin/env python3
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
'''

WORKER_SOURCE = r'''#!/bin/bash
set -euo pipefail
source "$HMP_CONDA_INIT"
conda activate "$HMP_ENV"
SID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$HMP_RUN/task_ids.txt")
[[ "$SID" =~ ^SRS[0-9]+$ ]] || { echo 'Invalid sample ID'; exit 1; }
if [ -s "$HMP_RUN/status/$SID.done" ]; then
  echo "Already complete: $SID"
  exit 0
fi
LOCK="$HMP_RUN/locks/$SID"
mkdir "$LOCK" || { echo "A task for $SID is already active"; exit 1; }
WORK=''
cleanup() {
  rc=$?
  set +e
  if [ "$rc" -ne 0 ]; then
    printf 'exit=%s\n' "$rc" > "$HMP_RUN/status/$SID.fail"
    if [ -n "$WORK" ] && [ -d "$WORK" ]; then
      [ ! -s "$WORK/report.tsv" ] || cp -n "$WORK/report.tsv" "$HMP_RUN/reports/$SID.failed.report.tsv"
      [ ! -s "$WORK/bracken.tsv" ] || cp -n "$WORK/bracken.tsv" "$HMP_RUN/reports/$SID.failed.bracken.tsv"
    fi
  fi
  cd "$HMP_RUN" 2>/dev/null || true
  if [ -n "$WORK" ] && [ -d "$WORK" ]; then
    python3 - "$HMP_RUN" "$WORK" "$SID" <<'PY'
import shutil, sys
from pathlib import Path
parent = (Path(sys.argv[1]) / 'work').resolve()
work = Path(sys.argv[2]).resolve()
if work.parent != parent or not work.name.startswith(sys.argv[3] + '_'):
    raise SystemExit('Refusing to remove unexpected temporary directory')
shutil.rmtree(work)
PY
  fi
  rmdir "$LOCK" 2>/dev/null || true
}
trap cleanup EXIT
WORK=$(mktemp -d "$HMP_RUN/work/${SID}_XXXXXX")
export TMPDIR="$WORK"
cd "$WORK"
echo "Sample: $SID; unpaired read mode; Bracken -r 100 -l S -t 10"
krakenuniq --version
bracken -v

python3 - "$HMP_BASE" "$HMP_RUN" "$SID" "$WORK" <<'PY' | sha256sum -c -
import csv, json, re, sys
from pathlib import Path
base, run, sid, work = map(Path, sys.argv[1:])
sid = str(sid)
with (run / 'manifest.tsv').open(newline='') as f:
    matches = [r for r in csv.DictReader(f, delimiter='\t') if r['sample'] == sid]
if len(matches) != 1 or matches[0]['depth'] != 'ge10000':
    raise SystemExit('Missing or ineligible sample')
row = matches[0]
audit = json.loads((base / 'audit/results/stool' / sid / 'AUDIT_COMPLETE.json').read_text())
if audit.get('site') != 'stool' or audit.get('SRS') != sid or audit.get('success') is not True:
    raise SystemExit('Audit identity/success mismatch')
root = (base / 'human/results/stool').resolve()
inputs = []
total = 0
for kind, col in [('R1_nohuman', 'R1'), ('R2_nohuman', 'R2'), ('singletons_nohuman', 'singleton')]:
    rec = audit['human'][kind]
    path = Path(rec['footprint']['path']).resolve()
    n = int(rec['reads'])
    if (n < 0 or n != int(row[col]) or path.name != f'{sid}_{kind}.fastq.gz'
        or path.parent not in (root, (root / sid).resolve())
        or rec.get('gzip_crc_valid') is not True
        or not re.fullmatch(r'[0-9a-f]{64}', rec.get('sha256', ''))):
        raise SystemExit(f'Invalid audit for {path}')
    print(f"{rec['sha256']}  {path}")
    total += n
    if n > 0:
        inputs.append(str(path))
if not inputs or total < 10000 or total != int(row['total_reads']):
    raise SystemExit('Read totals do not match the audited manifest')
(work / 'inputs.txt').write_text('\n'.join(inputs) + '\n')
(work / 'expected.txt').write_text(str(total) + '\n')
PY
mapfile -t INPUTS < "$WORK/inputs.txt"
[ "${#INPUTS[@]}" -gt 0 ] || { echo 'No input reads'; exit 1; }
EXPECTED=$(cat "$WORK/expected.txt")

python3 - "$WORK/inputs.txt" "$WORK/read_lengths.tsv" <<'PY'
from collections import Counter
from itertools import islice
from pathlib import Path
import gzip, sys
with open(sys.argv[2], 'x') as out:
    out.write('file\tsampled_reads\tlength_counts_top5\n')
    for filename in Path(sys.argv[1]).read_text().splitlines():
        with gzip.open(filename, 'rt') as f:
            counts = Counter(len(s.rstrip('\r\n')) for s in islice(f, 1, 4000, 4))
        out.write(f"{Path(filename).name}\t{sum(counts.values())}\t{counts.most_common(5)}\n")
PY

krakenuniq --db "$HMP_DB" --threads "$SLURM_CPUS_PER_TASK" --preload-size 48G \
  --report-file "$WORK/report.tsv" --output off "${INPUTS[@]}"
python3 - "$WORK/report.tsv" "$EXPECTED" <<'PY'
import csv, sys
vals = {}
with open(sys.argv[1], newline='') as f:
    for row in csv.reader(f, delimiter='\t'):
        if len(row) >= 9 and row[-3] in ('0', '1'):
            if row[-3] in vals:
                raise SystemExit('Duplicate root/unclassified rows')
            vals[row[-3]] = int(row[1])
if '1' not in vals or vals.get('0', 0) + vals['1'] != int(sys.argv[2]):
    raise SystemExit(f'Input/read count mismatch: {vals}; expected {sys.argv[2]}')
print('Read count verified:', vals.get('0', 0) + vals['1'])
PY
bracken -d "$HMP_DB" -i "$WORK/report.tsv" -o "$WORK/bracken.tsv" -r 100 -l S -t 10
test -s "$WORK/bracken.tsv"
python3 - "$WORK/bracken.tsv" <<'PY'
import csv, sys
with open(sys.argv[1], newline='') as f:
    r = csv.DictReader(f, delimiter='\t')
    if not r.fieldnames or not {'name','taxonomy_id','taxonomy_lvl','new_est_reads'} <= set(r.fieldnames):
        raise SystemExit('Unexpected Bracken output header')
    seen, total = set(), 0
    for row in r:
        tid = row['taxonomy_id']
        if not tid.isdigit() or tid in seen or row['taxonomy_lvl'] != 'S':
            raise SystemExit('Unexpected or duplicate species')
        seen.add(tid)
        n = int(row['new_est_reads'])
        if n < 0:
            raise SystemExit('Negative estimated read count')
        if tid not in ('9606','32630'):
            total += n
if total <= 0:
    raise SystemExit('No nonhost species reads')
print('Retained estimated species reads:', total)
PY
for name in report bracken read_lengths; do
  test ! -e "$HMP_RUN/reports/$SID.$name.tsv"
  mv "$WORK/$name.tsv" "$HMP_RUN/reports/$SID.$name.tsv"
done
rm -f -- "$HMP_RUN/status/$SID.fail"
printf 'KrakenUniq and Bracken passed; expected reads %s\n' "$EXPECTED" > "$HMP_RUN/status/$SID.done"
echo "Completed: $SID"
'''

MERGE_JOB_SOURCE = r'''#!/bin/bash
set -euo pipefail
python3 - "$HMP_RUN" "$HMP_PILOT_DIR" "$HMP_PILOT_SAMPLE" <<'PY'
import csv, os, shutil, sys
from pathlib import Path
run, pilot, sid = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]
if not (pilot / f'{sid}.done').is_file() or not (pilot / f'{sid}.bracken.tsv').is_file():
    raise SystemExit('Pilot result is missing')
with (run / 'manifest.tsv').open(newline='') as f:
    rows = list(csv.DictReader(f, delimiter='\t'))
done, excluded, problems = 0, 0, []
tmp = run / 'tables' / f'.sample_status.{os.getpid()}.tmp'
with tmp.open('x', newline='') as f:
    out = csv.writer(f, delimiter='\t', lineterminator='\n')
    out.writerow(['sample','depth','status'])
    for row in rows:
        sample, depth = row['sample'], row['depth']
        if depth == 'lt100':
            status = 'excluded'; excluded += 1
        elif sample == sid:
            status = 'done'; done += 1
        elif (run/'status'/f'{sample}.done').is_file() and (run/'reports'/f'{sample}.bracken.tsv').is_file():
            status = 'done'; done += 1
        else:
            status = 'fail' if (run/'status'/f'{sample}.fail').exists() else 'missing'
            problems.append(sample)
        out.writerow([sample, depth, status])
os.replace(tmp, run / 'tables/sample_status.tsv')
print(f'Profiles: {done}; excluded: {excluded}; failed/missing: {len(problems)}')
if done != 466 or excluded != 4 or problems:
    raise SystemExit('Incomplete run; see tables/sample_status.tsv')
for suffix, target in [('.bracken.tsv', '.bracken.tsv'), ('.krakenuniq.report.tsv', '.report.tsv')]:
    src, dst = pilot/f'{sid}{suffix}', run/'reports'/f'{sid}{target}'
    if not src.is_file():
        raise SystemExit(f'Missing pilot output: {src}')
    shutil.copyfile(src, dst)
PY
python3 "$HMP_RUN/scripts/merge.py" \
  --sample-status "$HMP_RUN/tables/sample_status.tsv" \
  --reports-dir "$HMP_RUN/reports" --output-dir "$HMP_RUN/tables" --expected-samples 466
'''


def run_command(args, env=None):
    return subprocess.run(args, env=env, text=True, capture_output=True, check=False)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--max-parallel', type=int, default=200)
    args = p.parse_args()
    if not 1 <= args.max_parallel <= 200:
        raise SystemExit('max-parallel must be between 1 and 200')

    base = Path(os.environ.get('HMP_BASE', '/fs/scratch/PAS1695/projects/oromic/shotgun/HMP_MV4_20261008/codex_autopilot_20261009_v1'))
    db = Path(os.environ.get('HMP_DB', '/fs/ess/PAS1695/globus/exorien_dbs/krakenuniq_exotic'))
    env_dir = Path(os.environ.get('HMP_ENV', '/users/PAS1695/nishiwaki/.conda/envs/krakenuniq_env'))
    conda_init = Path(os.environ.get('HMP_CONDA_INIT', '/apps/python/3.12/etc/profile.d/conda.sh'))
    meta = base/'stool_species_now_v425__a9bx868'
    pilot_dir = Path(os.environ.get('HMP_PILOT_DIR', str(base/'stool_kraken_bracken_pilot_SRS057049_KLD7yg')))
    pilot_id = '8301213'
    pilot_sample = 'SRS057049'
    for path in [db/'database.kdb',db/'database100mers.kmer_distrib',env_dir/'bin/krakenuniq',env_dir/'bin/bracken',conda_init,meta/'manifest.tsv',meta/'tables/sample_status.tsv']:
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f'Missing required file: {path}')
    with (meta/'manifest.tsv').open(newline='') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    with (meta/'tables/sample_status.tsv').open(newline='') as f:
        statuses = list(csv.DictReader(f, delimiter='\t'))
    if len(rows) != 470 or len(statuses) != 470:
        raise SystemExit('Unexpected manifest or status count')
    ids = [r['sample'] for r in rows]
    if len(set(ids)) != 470 or any(not re.fullmatch(r'SRS[0-9]+', s) for s in ids):
        raise SystemExit('Duplicate or invalid sample ID')
    status_map = {r['sample']: r['status'] for r in statuses}
    if len(status_map) != 470 or set(status_map) != set(ids):
        raise SystemExit('Sample status does not match manifest')
    selected = [r['sample'] for r in rows if r['depth'] == 'ge10000']
    excluded = [r['sample'] for r in rows if r['depth'] == 'lt100']
    if (len(selected) != 466 or len(excluded) != 4 or pilot_sample not in selected
        or any(status_map[s] != 'done' for s in selected)
        or any(status_map[s] != 'excluded' for s in excluded)):
        raise SystemExit('Unexpected sample selection or prior MetaPhlAn status')

    queued = run_command(['squeue','-h','-j',pilot_id,'-o','%T'])
    if queued.returncode != 0 and 'Invalid job id' not in queued.stderr:
        raise SystemExit('Unable to check pilot state: ' + queued.stderr.strip())
    live = queued.stdout.strip()
    if live:
        if live not in ('PENDING','RUNNING','COMPLETING','CONFIGURING'):
            raise SystemExit(f'Pilot not active in a usable state: {live}')
        dependency = f'afterok:{pilot_id}'
        print(f'Pilot {pilot_id}: {live}; the array will wait for success')
    else:
        accounting = run_command(['sacct','-n','-X','-j',pilot_id,'-P','--format=JobIDRaw,State,ExitCode'])
        if accounting.returncode != 0:
            raise SystemExit('Unable to check pilot completion: ' + accounting.stderr.strip())
        rec = next((line for line in accounting.stdout.splitlines() if line.split('|')[0].strip() == pilot_id), '')
        if [x.strip() for x in rec.split('|')] != [pilot_id, 'COMPLETED', '0:0']:
            raise SystemExit(f'Pilot not confirmed successful: {rec or "no accounting record"}')
        for path in (pilot_dir/f'{pilot_sample}.done', pilot_dir/f'{pilot_sample}.bracken.tsv', pilot_dir/f'{pilot_sample}.krakenuniq.report.tsv'):
            if not path.is_file() or path.stat().st_size == 0:
                raise SystemExit(f'Pilot result missing: {path}')
        dependency = None
        print(f'Pilot {pilot_id}: COMPLETED')

    lock = base/'stool_kraken_bracken_full_submission.lock'
    if lock.exists():
        raise SystemExit(f'Existing submission record: {lock}; inspect it before submitting again')
    with lock.open('x') as f:
        f.write('Preparing submission\n')
    run = None
    array_id = None
    submitted_unknown = False
    try:
        run = Path(tempfile.mkdtemp(prefix='stool_kraken_bracken_full_', dir=base))
        for name in ('scripts','logs','reports','status','locks','work','tables'):
            (run/name).mkdir()
        shutil.copyfile(meta/'manifest.tsv', run/'manifest.tsv')
        tasks = [sid for sid in selected if sid != pilot_sample]
        if len(tasks) != 465:
            raise RuntimeError('Unexpected array task count')
        (run/'task_ids.txt').write_text('\n'.join(tasks)+'\n')
        (run/'scripts/worker.sbatch').write_text(WORKER_SOURCE)
        (run/'scripts/merge.sbatch').write_text(MERGE_JOB_SOURCE)
        (run/'scripts/merge.py').write_text(MERGE_SOURCE)
        for path in (run/'scripts/worker.sbatch', run/'scripts/merge.sbatch'):
            chk = run_command(['bash','-n',str(path)])
            if chk.returncode != 0:
                raise RuntimeError(f'Script syntax check failed: {chk.stderr.strip()}')
        env = dict(os.environ, HMP_BASE=str(base), HMP_DB=str(db), HMP_ENV=str(env_dir),
                   HMP_CONDA_INIT=str(conda_init), HMP_RUN=str(run), HMP_PILOT_DIR=str(pilot_dir),
                   HMP_PILOT_SAMPLE=pilot_sample)
        common = ['sbatch','--parsable','--account=PAS1695','--partition=nextgen','--nodes=1']
        argv = common + ['--cpus-per-task=8','--mem=128G','--time=72:00:00',
                         f'--array=1-465%{args.max_parallel}','--job-name=hmp_gut_kb',
                         f'--output={run}/logs/profile_%A_%a.log']
        if dependency:
            argv += [f'--dependency={dependency}','--kill-on-invalid-dep=yes']
        argv += [str(run/'scripts/worker.sbatch')]
        result = run_command(argv, env=env)
        if result.returncode:
            raise RuntimeError('Array submission failed: '+result.stderr.strip())
        candidate = result.stdout.strip().split(';',1)[0]
        if not candidate.isdigit():
            submitted_unknown = True
            raise RuntimeError('Unknown array submission response: '+result.stdout.strip())
        array_id = candidate
        lock.write_text(f'Output directory: {run}\nArray job: {array_id}\nMerge: not yet submitted\n')
        print('Array job:', array_id, flush=True)
        argv = common + ['--cpus-per-task=2','--mem=16G','--time=04:00:00',
                         f'--dependency=afterany:{array_id}','--kill-on-invalid-dep=yes',
                         '--job-name=hmp_gut_kb_merge',f'--output={run}/logs/merge_%j.log',
                         str(run/'scripts/merge.sbatch')]
        result = run_command(argv, env=env)
        if result.returncode:
            raise RuntimeError('Merge submission failed: '+result.stderr.strip())
        merge_id = result.stdout.strip().split(';',1)[0]
        if not merge_id.isdigit():
            raise RuntimeError('Unknown merge submission response: '+result.stdout.strip())
        lock.write_text(f'Output directory: {run}\nArray job: {array_id}\nMerge job: {merge_id}\n')
        print('Merge job:', merge_id)
        print('Remaining samples:', len(tasks), '; pilot:', pilot_sample)
        print('Maximum simultaneous array tasks:', args.max_parallel)
        print('Output directory:', run)
    except Exception as exc:
        if array_id is None and not submitted_unknown:
            lock.unlink(missing_ok=True)
        else:
            with lock.open('a') as f:
                f.write(f'Submission incomplete: {exc}\n')
        raise SystemExit(f'ERROR: {exc}; inspect {lock} before any resubmission')


if __name__ == '__main__':
    main()
