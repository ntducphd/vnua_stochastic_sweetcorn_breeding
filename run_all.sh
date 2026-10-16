#!/usr/bin/env bash
# run_all.sh -- regenerate every table, figure and quoted number of the paper.
#
#   bash run_all.sh            from the deposited outputs (about one minute): pools the raw
#                              records, runs every analysis, draws the five figures, writes
#                              out/numbers.json, checks every quoted number and compares the
#                              regenerated files with checksums/outputs.sha256.
#   bash run_all.sh --full     runs the engine tests and reruns the simulations with sim/ from the
#                              deposited founder objects into out_full/ (about 13 CPU-hours; JOBS=4
#                              runs four jobs at a time) and compares every per-cycle record with
#                              the deposit.
#
# AlphaSimR 2.1.0 is required (see README); ALLOW_VERSION=1 skips the check.
set -euo pipefail
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1
MODE="${1:-outputs}"
JOBS="${JOBS:-4}"

ver=$(Rscript -e 'cat(as.character(packageVersion("AlphaSimR")))' 2>/dev/null)
if [ "$ver" != "2.1.0" ] && [ "${ALLOW_VERSION:-0}" != "1" ]; then
  echo "AlphaSimR $ver found; 2.1.0 is required (see README). Set ALLOW_VERSION=1 to override." >&2
  exit 1
fi

if [ "$MODE" = "--full" ]; then
  [ -e out_full ] && { echo "out_full/ exists; remove it first" >&2; exit 1; }
  echo "== engine tests (sim/tests/run_tests.R) =="
  Rscript sim/tests/run_tests.R
  echo "== Experiment 1 (the saved programme-layer founder object is reused; runMacs calls differ) =="
  mkdir -p out_full/exp1 && cp out/exp1/founder_20260925.rds* out_full/exp1/
  Rscript sim/run_exp1.R all --out out_full/exp1
  echo "== Experiment 2 (new founders, so Table 3 recurs only in distribution) =="
  mkdir -p out_full/exp2 out_full/application
  printf '%s\n' 0 1 2 3 | xargs -P "$JOBS" -I{} sh -c \
    'Rscript sim/run_exp2.R --shard {} --out out_full/exp2 > out_full/exp2_shard{}.log 2>&1'
  echo "== application, seeds 4001 to 4010, from the saved founder objects in out/founders/ =="
  seq 4001 4010 | xargs -P "$JOBS" -I{} sh -c \
    'Rscript sim/run_application.R --seed {} --out out_full/application --founders out/founders --pipeline-layout > out_full/app_{}.log 2>&1'
  python scripts/compare_rerun.py out_full/application
  exit 0
fi

echo "== pool the application records =="
python R/01_pool_application.py
echo "== Experiment 2: pool shards, variance-reduction table =="
python R/exp2_pool.py
python R/exp2_summary.py
echo "== application analyses =="
Rscript R/02_analyse_application.R
Rscript R/03_founder_scale.R
Rscript R/04_analyse_matched_controls.R
python R/05_bootstrap.py
echo "== inference -> out/numbers.json =="
python R/06_inference.py
echo "== figures =="
Rscript R/08_figures_3_5.R
Rscript R/07_figures_1_2.R
echo "== rerun of replicate 4001 =="
python scripts/compare_rerun.py out/rerun_4001
echo "== every quoted number =="
python scripts/verify_numbers.py | tail -1
echo "== regenerated files against checksums/outputs.sha256 =="
python - <<'PY'
import hashlib
bad = []
for line in open("checksums/outputs.sha256", encoding="utf-8"):
    want, path = line.split()
    path = path.lstrip("*")
    b = open(path, "rb").read()
    lf = b.replace(b"\r\n", b"\n")
    forms = (b,) if path.endswith(".png") else (lf, lf.replace(b"\n", b"\r\n"))
    if want not in {hashlib.sha256(x).hexdigest() for x in forms}:
        bad.append(path)
print("all regenerated files identical to the deposit" if not bad else "DIFFER: " + ", ".join(bad))
PY
python -c "import sys, numpy, scipy; print('this run: Python', sys.version.split()[0], '| NumPy', numpy.__version__, '| SciPy', scipy.__version__)"
