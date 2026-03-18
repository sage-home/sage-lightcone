#!/bin/bash

# Comprehensive validation tests for cli_lightcone argument checking.
# Each test runs cli_lightcone with a specific invalid (or valid) argument
# combination and checks that the expected error message appears on stderr,
# or that a valid combination is accepted without an error message.
#
# Exit codes from cli_lightcone are not yet reliable (validation errors
# currently exit 0), so all assertions are based on stderr content.

SCRIPT="${BASH_SOURCE[0]}"
[ -z "$SCRIPT" ] && SCRIPT="$0"
MY_SCRIPTS_DIRECTORY=$(cd "$(dirname "$SCRIPT")" && pwd)
MY_ROOT=$(cd "${MY_SCRIPTS_DIRECTORY}/.." && pwd)

CLI="${MY_ROOT}/bin/cli_lightcone"
DUMMY="tests/dummy.h5"   # does not need to exist for pre-file validation checks

PASS=0
FAIL=0

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# expect_error DESCRIPTION EXPECTED_PATTERN [cli args...]
# Passes if EXPECTED_PATTERN appears in stderr output.
expect_error() {
    local desc="$1"
    local pattern="$2"
    shift 2
    local output
    output=$("$CLI" "$@" 2>&1)
    if echo "$output" | grep -F -e "$pattern" -q; then
        echo "PASS: $desc"
        PASS=$((PASS + 1))
    else
        echo "FAIL: $desc"
        echo "      expected pattern: $pattern"
        echo "      actual output:    $output"
        FAIL=$((FAIL + 1))
    fi
}

# expect_no_error DESCRIPTION [cli args...]
# Passes if none of the known validation error strings appear in stderr output.
# (The run will still fail trying to open the non-existent dataset file — that's fine.)
expect_no_error() {
    local desc="$1"
    shift
    local output
    output=$("$CLI" "$@" 2>&1)
    if echo "$output" | grep -qE "Error: Invalid|Minimum DEC|Maximum DEC|Minimum RA|Maximum RA|Minimum redshift|Maximum redshift|outdir cannot be empty|outfile cannot be empty|filterfield cannot be empty|filtermin and|filtermax must be greater|must be numeric|dataset is required"; then
        echo "FAIL: $desc (unexpected validation error)"
        echo "      output: $output"
        FAIL=$((FAIL + 1))
    else
        echo "PASS: $desc"
        PASS=$((PASS + 1))
    fi
}

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------

cd "${MY_ROOT}" || exit 1

ORIGINAL_SCRIPTS_DIRECTORY=$MY_SCRIPTS_DIRECTORY
if [[ "$OSTYPE" == "darwin"* ]]; then
    source "${MY_ROOT}/setup_mac.sh"
elif [[ -f "/etc/modules" ]] || command -v module >/dev/null 2>&1; then
    source "${MY_ROOT}/setup.sh"
else
    source "${MY_ROOT}/setup.sh"
fi
MY_SCRIPTS_DIRECTORY=$ORIGINAL_SCRIPTS_DIRECTORY

if [ ! -f "$CLI" ]; then
    echo "ERROR: cli_lightcone not found at $CLI"
    echo "Run ./build_platform_aware.sh first."
    exit 1
fi

echo "==========================================="
echo "  cli_lightcone argument validation tests"
echo "==========================================="
echo ""

# ---------------------------------------------------------------------------
# --dataset
# ---------------------------------------------------------------------------
echo "--- dataset ---"
expect_error \
    "--dataset missing" \
    "--dataset is required" \
    --decmin 0 --decmax 1 --ramin 0 --ramax 1 --zmin 0 --zmax 1

# ---------------------------------------------------------------------------
# Declination
# ---------------------------------------------------------------------------
echo "--- declination ---"
expect_error \
    "decmin below -90" \
    "Minimum DEC cannot be less than -90" \
    --dataset $DUMMY --decmin -91 --decmax 0 --ramin 0 --ramax 1 --zmin 0 --zmax 1

expect_error \
    "decmin above 90" \
    "Minimum DEC cannot be less than -90 or greater than 90" \
    --dataset $DUMMY --decmin 100 --decmax 0 --ramin 0 --ramax 1 --zmin 0 --zmax 1

expect_error \
    "decmax below -90" \
    "Maximum DEC cannot be less than -90" \
    --dataset $DUMMY --decmin -90 --decmax -91 --ramin 0 --ramax 1 --zmin 0 --zmax 1

expect_error \
    "decmax above 90" \
    "Maximum DEC cannot be less than -90 or greater than 90" \
    --dataset $DUMMY --decmin 0 --decmax 100 --ramin 0 --ramax 1 --zmin 0 --zmax 1

expect_error \
    "decmin == decmax" \
    "Minimum DEC must be less than Maximum DEC" \
    --dataset $DUMMY --decmin 10 --decmax 10 --ramin 0 --ramax 1 --zmin 0 --zmax 1

expect_error \
    "decmin > decmax" \
    "Minimum DEC must be less than Maximum DEC" \
    --dataset $DUMMY --decmin 10 --decmax 0 --ramin 0 --ramax 1 --zmin 0 --zmax 1

# valid boundary values
expect_no_error \
    "dec at valid boundaries (-90 to 90)" \
    --dataset $DUMMY --decmin -90 --decmax 90 --ramin 0 --ramax 1 --zmin 0 --zmax 1

# ---------------------------------------------------------------------------
# Right ascension
# ---------------------------------------------------------------------------
echo "--- right ascension ---"
expect_error \
    "ramin below 0" \
    "Minimum RA cannot be less than zero" \
    --dataset $DUMMY --decmin 0 --decmax 1 --ramin -10 --ramax 1 --zmin 0 --zmax 1

expect_error \
    "ramin above 360" \
    "Minimum RA cannot be less than zero or greater than 360" \
    --dataset $DUMMY --decmin 0 --decmax 1 --ramin 400 --ramax 1 --zmin 0 --zmax 1

expect_error \
    "ramax below 0" \
    "Maximum RA cannot be less than zero" \
    --dataset $DUMMY --decmin 0 --decmax 1 --ramin 0 --ramax -1 --zmin 0 --zmax 1

expect_error \
    "ramax above 360" \
    "Maximum RA cannot be less than zero or greater than 360" \
    --dataset $DUMMY --decmin 0 --decmax 1 --ramin 0 --ramax 400 --zmin 0 --zmax 1

expect_error \
    "ramin == ramax" \
    "Minimum RA must be less than Maximum RA" \
    --dataset $DUMMY --decmin 0 --decmax 1 --ramin 10 --ramax 10 --zmin 0 --zmax 1

expect_error \
    "ramin > ramax" \
    "Minimum RA must be less than Maximum RA" \
    --dataset $DUMMY --decmin 0 --decmax 1 --ramin 10 --ramax 0 --zmin 0 --zmax 1

# valid boundary values
expect_no_error \
    "ra at valid boundaries (0 to 360)" \
    --dataset $DUMMY --decmin 0 --decmax 1 --ramin 0 --ramax 360 --zmin 0 --zmax 1

# ---------------------------------------------------------------------------
# Redshift
# ---------------------------------------------------------------------------
echo "--- redshift ---"
expect_error \
    "zmin negative" \
    "Minimum redshift must be greater than or equal to zero" \
    --dataset $DUMMY --decmin 0 --decmax 1 --ramin 0 --ramax 1 --zmin -1 --zmax 1

expect_error \
    "zmax == zmin" \
    "Minimum redshift must be less than Maximum redshift" \
    --dataset $DUMMY --decmin 0 --decmax 1 --ramin 0 --ramax 1 --zmin 1 --zmax 1

expect_error \
    "zmax < zmin" \
    "Minimum redshift must be less than Maximum redshift" \
    --dataset $DUMMY --decmin 0 --decmax 1 --ramin 0 --ramax 1 --zmin 3 --zmax 2

# valid boundary: zmin == 0
expect_no_error \
    "zmin == 0 is valid" \
    --dataset $DUMMY --decmin 0 --decmax 1 --ramin 0 --ramax 1 --zmin 0 --zmax 1

# ---------------------------------------------------------------------------
# Output directory
# ---------------------------------------------------------------------------
echo "--- output directory ---"
expect_error \
    "outdir empty string" \
    "--outdir cannot be empty" \
    --dataset $DUMMY --decmin 0 --decmax 1 --ramin 0 --ramax 1 --zmin 0 --zmax 1 \
    --outdir ""

# ---------------------------------------------------------------------------
# Output file
# ---------------------------------------------------------------------------
echo "--- output file ---"
expect_error \
    "outfile empty string" \
    "--outfile cannot be empty" \
    --dataset $DUMMY --decmin 0 --decmax 1 --ramin 0 --ramax 1 --zmin 0 --zmax 1 \
    --outfile ""

# ---------------------------------------------------------------------------
# Filter field
# ---------------------------------------------------------------------------
echo "--- filter ---"
expect_error \
    "filterfield empty string" \
    "--filterfield cannot be empty" \
    --dataset $DUMMY --decmin 0 --decmax 1 --ramin 0 --ramax 1 --zmin 0 --zmax 1 \
    --filterfield ""

expect_error \
    "filtermin and filtermax both empty" \
    "--filtermin and --filtermax cannot both be empty" \
    --dataset $DUMMY --decmin 0 --decmax 1 --ramin 0 --ramax 1 --zmin 0 --zmax 1 \
    --filtermin "" --filtermax ""

expect_error \
    "filtermin == filtermax" \
    "--filtermin and --filtermax cannot be the same" \
    --dataset $DUMMY --decmin 0 --decmax 1 --ramin 0 --ramax 1 --zmin 0 --zmax 1 \
    --filtermin 0 --filtermax 0

expect_error \
    "filtermax <= filtermin (filtermax < filtermin)" \
    "--filtermax must be greater than --filtermin" \
    --dataset $DUMMY --decmin 0 --decmax 1 --ramin 0 --ramax 1 --zmin 0 --zmax 1 \
    --filtermin 10 --filtermax 5

expect_error \
    "filtermax <= filtermin (filtermax == filtermin numeric)" \
    "--filtermin and --filtermax cannot be the same" \
    --dataset $DUMMY --decmin 0 --decmax 1 --ramin 0 --ramax 1 --zmin 0 --zmax 1 \
    --filtermin 5.0 --filtermax 5.0

expect_error \
    "filtermin non-numeric" \
    "--filtermin and --filtermax must be numeric values" \
    --dataset $DUMMY --decmin 0 --decmax 1 --ramin 0 --ramax 1 --zmin 0 --zmax 1 \
    --filtermin abc --filtermax 10

expect_error \
    "filtermax non-numeric" \
    "--filtermin and --filtermax must be numeric values" \
    --dataset $DUMMY --decmin 0 --decmax 1 --ramin 0 --ramax 1 --zmin 0 --zmax 1 \
    --filtermin 0 --filtermax abc

# ---------------------------------------------------------------------------
# --version
# ---------------------------------------------------------------------------
echo "--- version ---"
expect_no_error \
    "--version prints version without validation error" \
    --version

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo ""
echo "==========================================="
TOTAL=$((PASS + FAIL))
echo "Results: ${PASS}/${TOTAL} passed"
if [ $FAIL -gt 0 ]; then
    echo "  $FAIL test(s) FAILED"
    echo "==========================================="
    exit 1
else
    echo "All tests passed."
    echo "==========================================="
    exit 0
fi
