#!/bin/bash
(
set -euo pipefail
HMP_BASE="${HMP_BASE:-/fs/scratch/PAS1695/projects/oromic/shotgun/HMP_MV4_20261008/codex_autopilot_20261009_v1}"
HMP_DB="${HMP_DB:-/fs/ess/PAS1695/globus/exorien_dbs/krakenuniq_exotic}"
HMP_ENV="${HMP_ENV:-/users/PAS1695/nishiwaki/.conda/envs/krakenuniq_env}"
HMP_CONDA_INIT="${HMP_CONDA_INIT:-/apps/python/3.12/etc/profile.d/conda.sh}"
HMP_SAMPLE="SRS057049"
HMP_FASTQ="$HMP_BASE/human/results/stool/$HMP_SAMPLE"
HMP_META="$HMP_BASE/stool_species_now_v425__a9bx868"

test -s "$HMP_DB/database.kdb"
test -s "$HMP_DB/database100mers.kmer_distrib"
test -x "$HMP_ENV/bin/krakenuniq"
test -x "$HMP_ENV/bin/bracken"
test -s "$HMP_META/status/$HMP_SAMPLE.done"
for kind in R1_nohuman R2_nohuman singletons_nohuman; do
  test -s "$HMP_FASTQ/${HMP_SAMPLE}_${kind}.fastq.gz"
done
HMP_EXPECTED=$(python3 - "$HMP_META/manifest.tsv" "$HMP_SAMPLE" <<'PY'
import csv, sys
with open(sys.argv[1], newline='') as f:
    rows = [r for r in csv.DictReader(f, delimiter='\t') if r['sample'] == sys.argv[2]]
if len(rows) != 1 or rows[0]['depth'] != 'ge10000':
    raise SystemExit('Selected sample is absent or not eligible')
n = int(rows[0]['total_reads'])
if n < 10000:
    raise SystemExit('Sample has too few reads')
print(n)
PY
)
HMP_RUN=$(mktemp -d "$HMP_BASE/stool_kraken_bracken_pilot_${HMP_SAMPLE}_XXXXXX")
export HMP_BASE HMP_DB HMP_ENV HMP_CONDA_INIT HMP_SAMPLE HMP_FASTQ HMP_EXPECTED HMP_RUN

sbatch --account=PAS1695 --partition=nextgen --nodes=1 \
  --cpus-per-task=8 --mem=128G --time=24:00:00 \
  --job-name=hmp_gut_kb_pilot --output="$HMP_RUN/pilot_%j.log" <<'SBATCH'
#!/bin/bash
set -euo pipefail
source "$HMP_CONDA_INIT"
conda activate "$HMP_ENV"
echo "Sample: $HMP_SAMPLE; expected individual reads: $HMP_EXPECTED"
echo "DB: $HMP_DB; Bracken read-length model: 100 bp; unpaired read mode"
krakenuniq --version
bracken -v

python3 - "$HMP_BASE" "$HMP_SAMPLE" "$HMP_FASTQ" <<'PY' | sha256sum -c -
import json, re, sys
from pathlib import Path
base, sid, directory = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3])
audit = json.loads((base / 'audit/results/stool' / sid / 'AUDIT_COMPLETE.json').read_text())
if audit.get('site') != 'stool' or audit.get('SRS') != sid or audit.get('success') is not True:
    raise SystemExit('Audit identity/success mismatch')
for kind in ('R1_nohuman', 'R2_nohuman', 'singletons_nohuman'):
    rec = audit['human'][kind]
    path = directory / f'{sid}_{kind}.fastq.gz'
    if (Path(rec['footprint']['path']).resolve() != path.resolve()
        or rec.get('gzip_crc_valid') is not True
        or not re.fullmatch(r'[0-9a-f]{64}', rec.get('sha256', ''))):
        raise SystemExit(f'Invalid audit for {path}')
    print(f"{rec['sha256']}  {path}")
PY

krakenuniq --db "$HMP_DB" --threads "$SLURM_CPUS_PER_TASK" --preload-size 48G \
  --report-file "$HMP_RUN/$HMP_SAMPLE.krakenuniq.report.tsv" --output off \
  "$HMP_FASTQ/${HMP_SAMPLE}_R1_nohuman.fastq.gz" \
  "$HMP_FASTQ/${HMP_SAMPLE}_R2_nohuman.fastq.gz" \
  "$HMP_FASTQ/${HMP_SAMPLE}_singletons_nohuman.fastq.gz"

python3 - "$HMP_RUN/$HMP_SAMPLE.krakenuniq.report.tsv" "$HMP_EXPECTED" <<'PY'
import csv, sys
vals = {}
with open(sys.argv[1], newline='') as f:
    for row in csv.reader(f, delimiter='\t'):
        if len(row) >= 9 and row[-3] in ('0', '1'):
            if row[-3] in vals:
                raise SystemExit('Duplicate root/unclassified rows in KrakenUniq report')
            vals[row[-3]] = int(row[1])
if '1' not in vals or vals.get('0', 0) + vals['1'] != int(sys.argv[2]):
    raise SystemExit(f'Input/read count mismatch: {vals}; expected {sys.argv[2]}')
print('KrakenUniq read count verified:', vals.get('0', 0) + vals['1'])
PY

bracken -d "$HMP_DB" -i "$HMP_RUN/$HMP_SAMPLE.krakenuniq.report.tsv" \
  -o "$HMP_RUN/$HMP_SAMPLE.bracken.tsv" -r 100 -l S -t 10
test -s "$HMP_RUN/$HMP_SAMPLE.bracken.tsv"

python3 - "$HMP_RUN/$HMP_SAMPLE.bracken.tsv" "$HMP_RUN/$HMP_SAMPLE.nonhost_relative.tsv" "$HMP_SAMPLE" <<'PY'
import csv, sys
source, target, sample = sys.argv[1:]
with open(source, newline='') as f:
    reader = csv.DictReader(f, delimiter='\t')
    required = {'taxonomy_id', 'name', 'taxonomy_lvl', 'new_est_reads'}
    if not reader.fieldnames or not required.issubset(reader.fieldnames):
        raise SystemExit('Unexpected Bracken output columns')
    rows = [r for r in reader if r['taxonomy_id'] not in {'9606', '32630'}]
ids = [r['taxonomy_id'] for r in rows]
if len(ids) != len(set(ids)) or any(r['taxonomy_lvl'] != 'S' for r in rows):
    raise SystemExit('Unexpected or duplicated species in Bracken output')
counts = [int(r['new_est_reads']) for r in rows]
if any(c < 0 for c in counts) or sum(counts) <= 0:
    raise SystemExit('No valid retained species counts')
total = sum(counts)
with open(target, 'x', newline='') as f:
    writer = csv.writer(f, delimiter='\t', lineterminator='\n')
    writer.writerow(['taxonomy_id', 'name', 'taxonomy_lvl', sample])
    for row, count in zip(rows, counts):
        writer.writerow([row['taxonomy_id'], row['name'], row['taxonomy_lvl'], format(count / total, '.12g')])
print('Bracken species:', len(rows), 'retained estimated reads:', total)
print('Pilot normalized table:', target)
PY
printf 'KrakenUniq and Bracken pilot completed\n' > "$HMP_RUN/$HMP_SAMPLE.done"
SBATCH

echo "Pilot sample: $HMP_SAMPLE ($HMP_EXPECTED reads)"
echo "Pilot output: $HMP_RUN"
)
