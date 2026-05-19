#!/usr/bin/env bash
# Benchmark knock-knock 0.8.2 on 4 simulation quadrants under GNU /usr/bin/time -v.
# Pacbio/HiFi uses the CLI; Nanopore uses the run_nanopore.py driver (CLI doesn't dispatch nanopore).
set -euo pipefail

BENCH_DIR="$(cd "$(dirname "$0")"/.. && pwd)"
TIMEBIN=/home/gaox36/miniconda3/envs/py313/bin/time
VENV=/home/gaox36/bioinfo_tools/knock-knock/venv

# shellcheck disable=SC1091
source "$VENV/bin/activate"
export PATH="/home/gaox36/miniconda3/envs/py313/bin:$PATH"

TIMING="$BENCH_DIR/timing_logs"
mkdir -p "$TIMING"

echo "=== Benchmark run started: $(date -Iseconds) ==="
echo "Bench dir:   $BENCH_DIR"
echo "knock-knock: $(which knock-knock)"
echo "time:        $TIMEBIN"
echo "All 4 quadrants, 60,000 reads each, adapter-trimmed inputs"
echo

run_one() {
    local label="$1" batch="$2" sample="$3" mode="$4"
    local log="$TIMING/time_${label}.txt"
    local so="$TIMING/stdout_${label}.txt"
    local se="$TIMING/stderr_${label}.txt"
    echo "--- $label: start $(date -Iseconds) (mode=$mode) ---"
    cd "$BENCH_DIR/.."  # base dir so knock-knock CLI path args work
    if [[ "$mode" == "pacbio" ]]; then
        "$TIMEBIN" -v -o "$log" \
            knock-knock process-sample run_adapter_trim_bench "$batch" "$sample" \
                --stages preprocess,align,categorize \
                > "$so" 2> "$se" || { echo "FAILED: $label (exit $?)"; return 1; }
    else
        "$TIMEBIN" -v -o "$log" \
            "$VENV/bin/python3" "$BENCH_DIR/run_nanopore.py" run_adapter_trim_bench "$batch" "$sample" \
                > "$so" 2> "$se" || { echo "FAILED: $label (exit $?)"; return 1; }
    fi
    echo "--- $label: done $(date -Iseconds) ---"
    grep -E 'Elapsed|Maximum resident|Percent of CPU' "$log" | sed 's/^/    /'
    echo
}

# HiFi
run_one hifi_trac hifi_trac simulation_hifi_trac pacbio
run_one hifi_trbc hifi_trbc simulation_hifi_trbc pacbio

# ONT (uses driver)
run_one ont_trac  ont_trac  simulation_ont_trac  nanopore
run_one ont_trbc  ont_trbc  simulation_ont_trbc  nanopore

echo "=== Benchmark run finished: $(date -Iseconds) ==="
