#!/usr/bin/env bash
# Resume the 4-quadrant benchmark by running only the ONT samples.
# Fixes the python3 → venv python interpreter bug that caused the first run to fail.
set -euo pipefail

BENCH_DIR="$(cd "$(dirname "$0")"/.. && pwd)"
TIMEBIN=/home/gaox36/miniconda3/envs/py313/bin/time
VENV=/home/gaox36/bioinfo_tools/knock-knock/venv

# shellcheck disable=SC1091
source "$VENV/bin/activate"
export PATH="/home/gaox36/miniconda3/envs/py313/bin:$PATH"

TIMING="$BENCH_DIR/timing_logs"
mkdir -p "$TIMING"

echo "=== ONT-only bench resume: $(date -Iseconds) ==="
echo "Python:      $VENV/bin/python3"
echo "time:        $TIMEBIN"
echo

run_ont() {
    local label="$1" batch="$2" sample="$3"
    local log="$TIMING/time_${label}.txt"
    local so="$TIMING/stdout_${label}.txt"
    local se="$TIMING/stderr_${label}.txt"
    echo "--- $label: start $(date -Iseconds) ---"
    cd "$BENCH_DIR/.."
    "$TIMEBIN" -v -o "$log" \
        "$VENV/bin/python3" "$BENCH_DIR/run_nanopore.py" run_adapter_trim_bench "$batch" "$sample" \
            > "$so" 2> "$se" || { echo "FAILED: $label (exit $?)"; return 1; }
    echo "--- $label: done $(date -Iseconds) ---"
    grep -E 'Elapsed|Maximum resident|Percent of CPU' "$log" | sed 's/^/    /'
    echo
}

run_ont ont_trac ont_trac simulation_ont_trac
run_ont ont_trbc ont_trbc simulation_ont_trbc

echo "=== ONT bench finished: $(date -Iseconds) ==="
