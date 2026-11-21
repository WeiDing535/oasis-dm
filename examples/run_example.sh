#!/usr/bin/env bash
# Run oasis_dm on a bundled case and score the result against the reference phases.
#   ./examples/run_example.sh [tag]      # tag defaults to 7e1d; also 6zpp, 3o2e,
#                                         # and the NCS cases 3i0t, 3i0t_hard, 7tj1
#                                         # (an operator file tests/data/<tag>_ncs.txt
#                                         # is picked up automatically)
# Set HIST_LIB=/path/to/hist.lib if you want histogram matching; without it the
# script falls back to --no-hist.
set -euo pipefail
cd "$(dirname "$0")/.."

PY=${PYTHON:-python3}
TAG=${1:-7e1d}
CASE="tests/data/${TAG}.mtz"

# Run from the source tree even if the package is not pip-installed.
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"

if [ ! -f "$CASE" ]; then
  echo "ERROR: $CASE not found. Bundled tags:" >&2
  ls tests/data/*.mtz | sed 's|tests/data/||;s|\.mtz||' >&2
  exit 1
fi

mkdir -p out
ARGS=(--mtzin "$CASE" --mtzout "out/${TAG}_v2.mtz"
      --stats "out/${TAG}_v2_stats.json" --tag "$TAG" --solc 0.5)

# NCS cases ship with an operator file next to the data
if [ -f "tests/data/${TAG}_ncs.txt" ]; then
  ARGS+=(--ncs "tests/data/${TAG}_ncs.txt")
fi

if [ -n "${HIST_LIB:-}" ]; then
  ARGS+=(--hist-lib "$HIST_LIB")
elif ! "$PY" -c "import sys; sys.path.insert(0,'src'); from oasis_dm.paths import resolve_hist_lib; resolve_hist_lib()" 2>/dev/null; then
  echo "note: no hist.lib found (set HIST_LIB or CLIBD); running with --no-hist"
  ARGS+=(--no-hist)
fi

echo "=== [1/2] oasis_dm density modification: $TAG ==="
"$PY" -m oasis_dm "${ARGS[@]}"

echo "=== [2/2] score against reference phases ==="
"$PY" tools/evaluate.py \
  --case "$CASE" --truth-col PHIC --fo-col FOBS --fom-col FOM --free-col FreeR \
  --methods "initial=$CASE:PHIB:FOM,v2=out/${TAG}_v2.mtz:DM_PHIDM:DM_FOMDM" \
  --out "out/${TAG}_eval"

echo "done: out/${TAG}_v2.mtz, out/${TAG}_eval_metrics.json"
