#!/bin/bash
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
