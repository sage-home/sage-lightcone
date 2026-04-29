#!/bin/bash
set -e

# Locate project root so we can run from any directory
_F="${BASH_SOURCE[0]:-$0}"
_SCRIPT_DIR="$(cd "$(dirname "$_F")" && pwd)"
_PROJECT_ROOT="$(cd "${_SCRIPT_DIR}/.." && pwd)"

# Send output to log file and terminal (log file at project root)
MPI_TASKS=3
LOG_FILE="${_PROJECT_ROOT}/validate_kdtree_output_mpi${MPI_TASKS}.txt"
exec > >(tee -i ${LOG_FILE}) 2>&1

cd "${_PROJECT_ROOT}" || exit 1

# Validation script for sage2kdtree refactoring with MPI input

# Force MPI mode
export USE_MPI=yes

# 1. Setup Environment
echo "Detecting platform..."
if [[ "$OSTYPE" == "darwin"* ]]; then
    echo "macOS detected - using Homebrew setup"
    source setup_mac.sh
elif [[ -f "/etc/modules" ]] || command -v module >/dev/null 2>&1; then
    echo "HPC/Linux environment detected - assuming modules are loaded or sourcing setup.sh"
    source setup.sh || true
fi

# Verify executables were compiled with MPI — this script requires it.
# nm works on both Linux (ELF) and macOS (Mach-O); macOS prefixes C symbols
# with '_' but MPI_Init still matches as a substring either way.
for _bin in "${_PROJECT_ROOT}/bin/sage2kdtree" "${_PROJECT_ROOT}/bin/cli_lightcone"; do
    if ! nm "$_bin" 2>/dev/null | grep -q 'MPI_Init'; then
        echo "ERROR: $(basename $_bin) was not built with MPI support."
        echo "       Rebuild with MPI enabled first:"
        echo "         USE_MPI=yes ./build_platform_aware.sh"
        exit 1
    fi
done

# 2. Define Inputs
TEST_DIR="tests/sage-model-tests"

# Check for benchmark output first (CI default), then standard output (local default)
INPUT_SAGE_DIR="${TEST_DIR}/output_sage_hdf5_mpi${MPI_TASKS}/millennium"
echo "Using SAGE output: ${INPUT_SAGE_DIR}"

PARAM_FILE="${TEST_DIR}/input/millennium_sage_hdf5.par"
ALIST_FILE="${TEST_DIR}/input/millennium/trees/millennium.a_list"

# 3. Validation Output Directory
VALIDATION_DIR="validation_output"
#rm -rf ${VALIDATION_DIR}
mkdir -p ${VALIDATION_DIR}

echo "Using SAGE input from: ${INPUT_SAGE_DIR}"

if [ ! -f "${INPUT_SAGE_DIR}/model_0.hdf5" ]; then
    echo "ERROR: SAGE output not found at ${INPUT_SAGE_DIR}/model_0.hdf5"
    echo "Please run 'tests/test_sage_hdf5_mpi.sh --np ${MPI_TASKS}' first to generate input data."
    exit 1
fi



# 4. Run New Code (sage2kdtree)
echo ""
echo "=========================================================="
echo "Running NEW sage2kdtree (the test)..."
echo "=========================================================="
start_new=$SECONDS
./bin/sage2kdtree \
    -s ${INPUT_SAGE_DIR} \
    -p ${PARAM_FILE} \
    -a ${ALIST_FILE} \
    -o ${VALIDATION_DIR}/kdtree_test_mpi${MPI_TASKS}.h5 \
    --ppc 1000 \
    -v 1
duration_new=$(( SECONDS - start_new ))
echo "NEW sage2kdtree took ${duration_new} seconds."

# Run cli_lightcone to generate lightcone_3d_SnapNum.png for visual verification
echo "========== PHASE 3: Running cli_lightcone =========="
mpirun -np $MPI_TASKS ./bin/cli_lightcone \
    --dataset ${VALIDATION_DIR}/kdtree_test_mpi${MPI_TASKS}.h5 \
    --decmin 0 --decmax 1 --ramin 0 --ramax 1 --zmin 0 --zmax 1 \
    --outdir ${VALIDATION_DIR} --outfile lightcone.h5
    #--outdir ${OUTPUTDIR} --outfields ra dec redshift_cosmological redshift_observed sfr --outfile $RAWNAME-lightcone.h5


# 5: Plot (not benchmarked - visualization only)
echo "========== PHASE 4: Plotting =========="
source .venv/bin/activate
# Field names are case-insensitive in plot_lightcone.py (SnapNum = snapnum)
python3 ./src/plot_lightcone.py ${VALIDATION_DIR}/lightcone.h5 SnapNum