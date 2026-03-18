#!/bin/bash
#
# profile_xctrace.sh
#
# Profiles sage2kdtree from both 'main' and 'before' branches using xctrace
# (Time Profiler) and the 'sample' command (plain-text call tree).
#
# Usage:
#   ./profile_xctrace.sh
#
# Output files (in ./profile_output/):
#   main_sage2kdtree.trace      - open in Instruments.app
#   before_sage2kdtree.trace    - open in Instruments.app
#   main_sample.txt             - plain-text call tree
#   before_sample.txt           - plain-text call tree

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAIN_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
BEFORE_ROOT="$(cd "${MAIN_ROOT}/../before" && pwd)"
PROFILE_DIR="${MAIN_ROOT}/tests/sage-model-tests/profile_output"

MAIN_S2K="${MAIN_ROOT}/bin/sage2kdtree"
BEFORE_S2K="${BEFORE_ROOT}/bin/sage2kdtree"
SAGE_BIN="${MAIN_ROOT}/bin/sage"

SAGE_OUTPUT="${PROFILE_DIR}/sage_output"
S2K_ARGS_COMMON="-s ${SAGE_OUTPUT} -p input/millennium_sage_hdf5.par -a input/millennium/trees/millennium.a_list --ppc 1000 -v 1 --noarrays"

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------

cd "${MAIN_ROOT}/tests/sage-model-tests"

# Source platform environment (sets HDF5/Boost/etc paths)
if [ -f "${MAIN_ROOT}/setup_mac.sh" ]; then
    source "${MAIN_ROOT}/setup_mac.sh" 2>/dev/null || true
fi

echo "Checking binaries..."
for bin in "${MAIN_S2K}" "${BEFORE_S2K}" "${SAGE_BIN}"; do
    if [ ! -f "${bin}" ]; then
        echo "ERROR: Missing binary: ${bin}"
        exit 1
    fi
done
echo "  main:   ${MAIN_S2K}"
echo "  before: ${BEFORE_S2K}"

mkdir -p "${PROFILE_DIR}"

# ---------------------------------------------------------------------------
# Run SAGE if needed
# ---------------------------------------------------------------------------

if [ ! -d "${SAGE_OUTPUT}" ]; then
    echo ""
    echo "========== Running SAGE (shared across all profiles) =========="
    # Prepare parameter files (reuse test script helpers)
    if [ ! -f input/millennium.par ]; then
        cp sage-model/input/millennium.par input/millennium.par 2>/dev/null || \
        cp "${MAIN_ROOT}/sage-model/input/millennium.par" input/millennium.par
    fi
    # Run SAGE
    "${SAGE_BIN}" input/millennium_sage_hdf5.par
    # Move output to profile_output/sage_output
    mv output/millennium "${SAGE_OUTPUT}"
    echo "SAGE complete. Output in: ${SAGE_OUTPUT}"
else
    echo "SAGE output already exists: ${SAGE_OUTPUT} (reusing)"
fi

# ---------------------------------------------------------------------------
# Profile helper
# ---------------------------------------------------------------------------

run_profile() {
    local label="$1"
    local binary="$2"
    local trace_file="${PROFILE_DIR}/${label}_sage2kdtree.trace"
    local sample_file="${PROFILE_DIR}/${label}_sample.txt"
    local output_kdtree="${PROFILE_DIR}/${label}_kdtree.h5"

    echo ""
    echo "========== Profiling: ${label} =========="
    echo "  Binary:  ${binary}"
    echo "  Trace:   ${trace_file}"
    echo "  Sample:  ${sample_file}"

    # Remove previous trace if it exists (xctrace refuses to overwrite)
    rm -rf "${trace_file}"

    # Run xctrace (Time Profiler)
    echo "  Running xctrace..."
    xctrace record \
        --template 'Time Profiler' \
        --output "${trace_file}" \
        --launch -- \
        "${binary}" \
        ${S2K_ARGS_COMMON} \
        -o "${output_kdtree}" \
        2>&1 | grep -v "^$" | head -10

    echo "  xctrace complete: ${trace_file}"

    # Export heavy-stack summary to text
    echo "  Exporting call tree summary..."
    xctrace export \
        --input "${trace_file}" \
        --xpath '//trace-toc/run[@number="1"]/data/table[@schema="time-profile"]' \
        2>/dev/null | python3 - << 'PYEOF'
import sys, xml.etree.ElementTree as ET
data = sys.stdin.read()
if not data.strip():
    print("(no data exported)")
    sys.exit(0)
try:
    root = ET.fromstring(data)
    print("exported ok, rows:", len(list(root)))
except Exception as e:
    print("XML parse error:", e)
    print(data[:500])
PYEOF
    echo "  (Full trace in ${trace_file} — open in Instruments.app for interactive view)"
}

run_sample() {
    local label="$1"
    local binary="$2"
    local sample_file="${PROFILE_DIR}/${label}_sample.txt"
    local output_kdtree="${PROFILE_DIR}/${label}_kdtree_sample.h5"

    echo ""
    echo "========== Sampling: ${label} =========="

    # Run sage2kdtree in background, sample it, wait
    rm -f "${output_kdtree}"
    "${binary}" \
        ${S2K_ARGS_COMMON} \
        -o "${output_kdtree}" \
        > /dev/null 2>&1 &
    local pid=$!
    echo "  PID ${pid} started, sampling for up to 60s..."
    sample "${pid}" 60 -file "${sample_file}" -wait 2>/dev/null || true
    wait "${pid}" 2>/dev/null || true
    echo "  Sample saved: ${sample_file}"
}

# ---------------------------------------------------------------------------
# Run profiles
# ---------------------------------------------------------------------------

echo ""
echo "NOTE: Using -v 1 (not -v 2) to reduce I/O noise in profiles."
echo ""

run_sample "main"   "${MAIN_S2K}"
run_sample "before" "${BEFORE_S2K}"
run_profile "main"   "${MAIN_S2K}"
run_profile "before" "${BEFORE_S2K}"

# ---------------------------------------------------------------------------
# Quick comparison: top functions from sample output
# ---------------------------------------------------------------------------

echo ""
echo "=================================================================="
echo "TOP FUNCTIONS BY SAMPLE COUNT"
echo "=================================================================="

for label in main before; do
    sample_file="${PROFILE_DIR}/${label}_sample.txt"
    if [ -f "${sample_file}" ]; then
        echo ""
        echo "--- ${label} ---"
        # Extract non-blank, non-header lines that look like call counts
        grep -E "^\s+[0-9]+ " "${sample_file}" 2>/dev/null | \
            sort -rn | head -30 || \
        grep -v "^Call graph\|^Binary\|^Date\|^OS X\|^Start\|^End\|^Sampling\|^Parent\|^Process\|^[[:space:]]*$" \
            "${sample_file}" | head -40
    fi
done

echo ""
echo "=================================================================="
echo "Trace files for Instruments.app:"
echo "  open ${PROFILE_DIR}/main_sage2kdtree.trace"
echo "  open ${PROFILE_DIR}/before_sage2kdtree.trace"
echo ""
echo "Or diff sample outputs:"
echo "  diff ${PROFILE_DIR}/main_sample.txt ${PROFILE_DIR}/before_sample.txt | head -100"
echo "=================================================================="
