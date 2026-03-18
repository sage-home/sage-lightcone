#!/bin/bash

# Comprehensive validation tests for sage2kdtree argument checking.
# Each test runs sage2kdtree with a specific invalid (or valid) argument
# combination and checks that the expected error message appears on stderr,
# or that a valid combination is accepted without an error message.

SCRIPT="${BASH_SOURCE[0]}"
[ -z "$SCRIPT" ] && SCRIPT="$0"
MY_SCRIPTS_DIRECTORY=$(cd "$(dirname "$SCRIPT")" && pwd)
MY_ROOT=$(cd "${MY_SCRIPTS_DIRECTORY}/.."; pwd)

CLI="${MY_ROOT}/bin/sage2kdtree"

PASS=0
FAIL=0

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# expect_error DESCRIPTION EXPECTED_PATTERN [cli args...]
# Passes if EXPECTED_PATTERN appears in combined stdout+stderr output.
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
# Passes if none of the known argument-error strings appear in output.
# (The run will still fail trying to open non-existent paths — that is fine.)
expect_no_error() {
    local desc="$1"
    shift
    local output
    output=$("$CLI" "$@" 2>&1)
    if echo "$output" | grep -qE "is required but missing|unrecognised option"; then
        echo "FAIL: $desc (unexpected argument error)"
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
    echo "ERROR: sage2kdtree not found at $CLI"
    echo "Run ./build_platform_aware.sh first."
    exit 1
fi

# Dummy paths — do not need to exist for argument-presence checks.
DUMMY_DIR="/tmp/sage_dummy_dir"
DUMMY_FILE="/tmp/sage_dummy.par"
DUMMY_ALIST="/tmp/sage_dummy_alist.txt"
DUMMY_OUTPUT="/tmp/sage_dummy_output.h5"

echo "==========================================="
echo "  sage2kdtree argument validation tests"
echo "==========================================="
echo ""

# ---------------------------------------------------------------------------
# Required arguments — each missing individually
# ---------------------------------------------------------------------------
echo "--- required arguments ---"

expect_error \
    "--sage missing" \
    "the option '--sage' is required but missing" \
    -p "$DUMMY_FILE" -a "$DUMMY_ALIST" -o "$DUMMY_OUTPUT"

expect_error \
    "--param missing" \
    "the option '--param' is required but missing" \
    -s "$DUMMY_DIR" -a "$DUMMY_ALIST" -o "$DUMMY_OUTPUT"

expect_error \
    "--alist missing" \
    "the option '--alist' is required but missing" \
    -s "$DUMMY_DIR" -p "$DUMMY_FILE" -o "$DUMMY_OUTPUT"

expect_error \
    "--output missing" \
    "the option '--output' is required but missing" \
    -s "$DUMMY_DIR" -p "$DUMMY_FILE" -a "$DUMMY_ALIST"

expect_error \
    "all required args missing" \
    "is required but missing" \

# ---------------------------------------------------------------------------
# Unrecognised options
# ---------------------------------------------------------------------------
echo "--- unrecognised options ---"

expect_error \
    "unrecognised long option" \
    "unrecognised option" \
    -s "$DUMMY_DIR" -p "$DUMMY_FILE" -a "$DUMMY_ALIST" -o "$DUMMY_OUTPUT" \
    --notanoption

expect_error \
    "unrecognised short option" \
    "unrecognised option" \
    -s "$DUMMY_DIR" -p "$DUMMY_FILE" -a "$DUMMY_ALIST" -o "$DUMMY_OUTPUT" \
    -Z

# ---------------------------------------------------------------------------
# Valid invocation (all required args present; paths need not exist —
# argument parsing succeeds, but file-not-found errors are acceptable)
# ---------------------------------------------------------------------------
echo "--- valid argument combinations ---"

expect_no_error \
    "all required args present (paths may not exist)" \
    -s "$DUMMY_DIR" -p "$DUMMY_FILE" -a "$DUMMY_ALIST" -o "$DUMMY_OUTPUT"

expect_no_error \
    "all required args with optional --ppc" \
    -s "$DUMMY_DIR" -p "$DUMMY_FILE" -a "$DUMMY_ALIST" -o "$DUMMY_OUTPUT" \
    --ppc 500

expect_no_error \
    "all required args with --centralgalaxies flag" \
    -s "$DUMMY_DIR" -p "$DUMMY_FILE" -a "$DUMMY_ALIST" -o "$DUMMY_OUTPUT" \
    --centralgalaxies

expect_no_error \
    "all required args with --noarrays flag" \
    -s "$DUMMY_DIR" -p "$DUMMY_FILE" -a "$DUMMY_ALIST" -o "$DUMMY_OUTPUT" \
    --noarrays

expect_no_error \
    "long-form required arg names" \
    --sage "$DUMMY_DIR" --param "$DUMMY_FILE" --alist "$DUMMY_ALIST" --output "$DUMMY_OUTPUT"

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
