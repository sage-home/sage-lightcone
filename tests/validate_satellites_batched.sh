#!/bin/bash
#
# validate_satellites_batched.sh
#
# Validates that --satellites_batched produces identical output to the default
# (unbatched) satellite read path in cli_lightcone.
#
# Runs cli_lightcone twice against the same KD-tree file:
#   Run 1: default (unbatched) satellite reads
#   Run 2: --satellites_batched
#
# Both outputs are sorted by (SnapNum, GalaxyIndex, ra, dec) before comparison
# so that any ordering differences introduced by the batched read path do not
# cause false failures.
#
# Usage:
#   ./tests/validate_satellites_batched.sh [--dataset <kdtree.h5>] [--rebuild]
#
# If --dataset is omitted the script tries to find a suitable KD-tree file
# automatically (searching validation_output/ and tests/sage-model-tests/).

# ---------------------------------------------------------------------------
# Timing helper
# ---------------------------------------------------------------------------

now_sec() { python3 -c "import time; print(f'{time.time():.3f}')"; }
elapsed_sec() { python3 -c "print(f'{$2 - $1:.1f}')"; }

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------

FORCE_REBUILD=0
DATASET=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --dataset)
            DATASET="$2"; shift 2 ;;
        --rebuild)
            FORCE_REBUILD=1; shift ;;
        --help|-h)
            echo "Usage: $0 [--dataset <kdtree.h5>] [--rebuild]"
            echo ""
            echo "  --dataset <file>  KD-tree HDF5 file to use as input."
            echo "                    If omitted, the script searches for one automatically."
            echo "  --rebuild         Force rebuild of executables before running."
            exit 0 ;;
        *)
            echo "Unknown option: $1"; exit 1 ;;
    esac
done

# ---------------------------------------------------------------------------
# Locate root and set up environment
# ---------------------------------------------------------------------------

MY_SCRIPT="${BASH_SOURCE[0]:-$0}"
MY_SCRIPTS_DIRECTORY=$(cd "$(dirname "$MY_SCRIPT")" && pwd)
MY_ROOT=$(cd "${MY_SCRIPTS_DIRECTORY}/.." && pwd)

ORIGINAL_SCRIPTS_DIRECTORY=$MY_SCRIPTS_DIRECTORY
if [[ "$OSTYPE" == "darwin"* ]]; then
    source "${MY_ROOT}/setup_mac.sh"
elif [[ -f "/etc/modules" ]] || command -v module >/dev/null 2>&1; then
    source "${MY_ROOT}/setup.sh"
else
    source "${MY_ROOT}/setup.sh"
fi
MY_SCRIPTS_DIRECTORY=$ORIGINAL_SCRIPTS_DIRECTORY

CLI="${MY_ROOT}/bin/cli_lightcone"

# ---------------------------------------------------------------------------
# Optionally rebuild
# ---------------------------------------------------------------------------

if [ $FORCE_REBUILD -eq 1 ] || [ ! -f "$CLI" ]; then
    echo "Building executables..."
    pushd "${MY_ROOT}" > /dev/null
    ./build_platform_aware.sh
    BUILD_STATUS=$?
    popd > /dev/null
    if [ $BUILD_STATUS -ne 0 ]; then
        echo "ERROR: Build failed"
        exit 1
    fi
fi

if [ ! -f "$CLI" ]; then
    echo "ERROR: cli_lightcone not found at $CLI"
    echo "Run ./build_platform_aware.sh first, or pass --rebuild."
    exit 1
fi

# ---------------------------------------------------------------------------
# Locate KD-tree dataset
# ---------------------------------------------------------------------------

if [ -z "$DATASET" ]; then
    # Try common locations
    CANDIDATES=(
        "${MY_ROOT}/validation_output/kdtree_test_mpi3.h5"
        "${MY_ROOT}/tests/sage-model-tests/output/myhdf5millennium-kdtree.h5"
        "${MY_ROOT}/tests/sage-model-tests/output_sage_hdf5_one_step/myhdf5millennium-kdtree-onestep.h5"
    )
    for c in "${CANDIDATES[@]}"; do
        if [ -f "$c" ]; then
            DATASET="$c"
            echo "Using dataset: $DATASET"
            break
        fi
    done
fi

if [ -z "$DATASET" ] || [ ! -f "$DATASET" ]; then
    echo "ERROR: No KD-tree HDF5 dataset found."
    echo "Run tests/validate_kdtree_output.sh first to produce one, or pass --dataset <file>."
    exit 1
fi

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

OUTDIR="${MY_ROOT}/tests/validate_satellites_batched_output"
mkdir -p "$OUTDIR"

OUT_UNBATCHED="${OUTDIR}/lightcone_unbatched.h5"
OUT_BATCHED="${OUTDIR}/lightcone_batched.h5"

LC_DECMIN=0; LC_DECMAX=1
LC_RAMIN=0;  LC_RAMAX=1
LC_ZMIN=0;   LC_ZMAX=1

COMMON_ARGS=(
    --dataset "$DATASET"
    --centralgalaxies
    --decmin $LC_DECMIN --decmax $LC_DECMAX
    --ramin  $LC_RAMIN  --ramax  $LC_RAMAX
    --zmin   $LC_ZMIN   --zmax   $LC_ZMAX
    --outdir "$OUTDIR"
)

# ---------------------------------------------------------------------------
# Run 1: default (unbatched)
# ---------------------------------------------------------------------------

echo ""
echo "========== Run 1: default (unbatched satellite reads) =========="
T1_START=$(now_sec)
"$CLI" "${COMMON_ARGS[@]}" --outfile lightcone_unbatched.h5
RUN1_STATUS=$?
T1_END=$(now_sec)
T1=$(elapsed_sec $T1_START $T1_END)

if [ $RUN1_STATUS -ne 0 ] || [ ! -f "$OUT_UNBATCHED" ]; then
    echo "ERROR: Run 1 (unbatched) failed (exit code $RUN1_STATUS)"
    exit 1
fi
echo "Run 1 complete in ${T1}s"

# ---------------------------------------------------------------------------
# Run 2: batched
# ---------------------------------------------------------------------------

echo ""
echo "========== Run 2: --satellites_batched =========="
T2_START=$(now_sec)
"$CLI" "${COMMON_ARGS[@]}" --satellites_batched --outfile lightcone_batched.h5
RUN2_STATUS=$?
T2_END=$(now_sec)
T2=$(elapsed_sec $T2_START $T2_END)

if [ $RUN2_STATUS -ne 0 ] || [ ! -f "$OUT_BATCHED" ]; then
    echo "ERROR: Run 2 (batched) failed (exit code $RUN2_STATUS)"
    exit 1
fi
echo "Run 2 complete in ${T2}s"

# ---------------------------------------------------------------------------
# Compare outputs (sorted)
# ---------------------------------------------------------------------------

echo ""
echo "========== Comparing outputs =========="

COMPARE_RESULT=$("${MY_ROOT}/.venv/bin/python3" - "$OUT_UNBATCHED" "$OUT_BATCHED" <<'PYEOF'
import sys
import h5py
import numpy as np

def load_and_sort(path):
    """Load all root-level datasets from a lightcone HDF5 file and return them
    sorted by (SnapNum, GalaxyIndex, ra, dec) to normalise row ordering."""
    with h5py.File(path, 'r') as f:
        # Collect all root-level scalar (1D) and array (2D) datasets
        scalar_fields = {}
        array_fields = {}
        for key in f.keys():
            obj = f[key]
            if not isinstance(obj, h5py.Dataset):
                continue
            if len(obj.shape) == 1:
                scalar_fields[key] = obj[:]
            elif len(obj.shape) == 2:
                array_fields[key] = obj[:]

    if not scalar_fields:
        return None, None, None, "No scalar datasets found"

    # Determine galaxy count from any 1D field
    n = next(iter(scalar_fields.values())).shape[0]

    # Build sort key: prefer (SnapNum, GalaxyIndex, ra, dec); fall back gracefully
    def find_key(fields, candidates):
        for c in candidates:
            for k in fields:
                if k.lower() == c.lower():
                    return fields[k]
        return None

    snap = find_key(scalar_fields, ['SnapNum', 'snapnum'])
    gidx = find_key(scalar_fields, ['GalaxyIndex', 'galaxyindex'])
    ra   = find_key(scalar_fields, ['ra', 'Ra', 'RA'])
    dec  = find_key(scalar_fields, ['dec', 'Dec', 'DEC'])

    sort_parts = [x for x in [snap, gidx, ra, dec] if x is not None]
    if sort_parts:
        # Build structured array for lexsort (lexsort sorts by last key first)
        sort_keys = tuple(reversed(sort_parts))
        order = np.lexsort(sort_keys)
    else:
        order = np.arange(n)

    sorted_scalars = {k: v[order] for k, v in scalar_fields.items()}
    sorted_arrays  = {k: v[order] for k, v in array_fields.items()}
    return sorted_scalars, sorted_arrays, n, None

file1, file2 = sys.argv[1], sys.argv[2]

scalars1, arrays1, n1, err1 = load_and_sort(file1)
scalars2, arrays2, n2, err2 = load_and_sort(file2)

if err1:
    print(f"FAIL: could not load {file1}: {err1}")
    sys.exit(1)
if err2:
    print(f"FAIL: could not load {file2}: {err2}")
    sys.exit(1)

failures = []

# Check galaxy count
if n1 != n2:
    failures.append(f"galaxy count mismatch: unbatched={n1}, batched={n2}")
else:
    print(f"  Galaxy count: {n1}  (match)")

# Check scalar fields
keys1 = set(scalars1.keys())
keys2 = set(scalars2.keys())
if keys1 != keys2:
    only1 = keys1 - keys2
    only2 = keys2 - keys1
    if only1:
        failures.append(f"fields only in unbatched: {sorted(only1)}")
    if only2:
        failures.append(f"fields only in batched:   {sorted(only2)}")

for key in sorted(keys1 & keys2):
    a = scalars1[key]
    b = scalars2[key]
    if a.shape != b.shape:
        failures.append(f"{key}: shape mismatch {a.shape} vs {b.shape}")
        continue
    if np.issubdtype(a.dtype, np.floating):
        if not np.allclose(a, b, rtol=1e-6, atol=1e-12, equal_nan=True):
            max_diff = np.max(np.abs(a - b))
            failures.append(f"{key}: float mismatch (max diff={max_diff:.3e})")
        else:
            print(f"  {key}: OK")
    else:
        if not np.array_equal(a, b):
            n_diff = int(np.sum(a != b))
            failures.append(f"{key}: integer mismatch ({n_diff}/{len(a)} elements differ)")
        else:
            print(f"  {key}: OK")

# Check array (2D) fields
akeys1 = set(arrays1.keys())
akeys2 = set(arrays2.keys())
for key in sorted(akeys1 & akeys2):
    a = arrays1[key]
    b = arrays2[key]
    if a.shape != b.shape:
        failures.append(f"{key} (2D): shape mismatch {a.shape} vs {b.shape}")
        continue
    if np.issubdtype(a.dtype, np.floating):
        if not np.allclose(a, b, rtol=1e-6, atol=1e-12, equal_nan=True):
            max_diff = float(np.max(np.abs(a - b)))
            failures.append(f"{key} (2D): float mismatch (max diff={max_diff:.3e})")
        else:
            print(f"  {key} (2D): OK")
    else:
        if not np.array_equal(a, b):
            n_diff = int(np.sum(a != b))
            failures.append(f"{key} (2D): integer mismatch ({n_diff} elements differ)")
        else:
            print(f"  {key} (2D): OK")

if failures:
    print("")
    for f in failures:
        print(f"  FAIL: {f}")
    sys.exit(1)
else:
    sys.exit(0)
PYEOF
)

COMPARE_STATUS=$?
echo "$COMPARE_RESULT"

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

echo ""
echo "======================================================================"
echo "  SATELLITES_BATCHED VALIDATION REPORT"
echo "======================================================================"
echo ""
echo "  Dataset:   $DATASET"
echo "  Lightcone: dec=[${LC_DECMIN},${LC_DECMAX}] ra=[${LC_RAMIN},${LC_RAMAX}] z=[${LC_ZMIN},${LC_ZMAX}]"
echo ""
printf "  %-30s  %8ss\n" "Run 1 (unbatched)" "${T1}"
printf "  %-30s  %8ss\n" "Run 2 (--satellites_batched)" "${T2}"
echo ""

if [ $COMPARE_STATUS -eq 0 ]; then
    echo "  Result: PASS — both outputs are identical (after sorting)"
    echo "======================================================================"
    exit 0
else
    echo "  Result: FAIL — outputs differ (see above)"
    echo "======================================================================"
    exit 1
fi
