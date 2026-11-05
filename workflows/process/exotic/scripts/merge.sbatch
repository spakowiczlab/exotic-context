#!/bin/bash
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
