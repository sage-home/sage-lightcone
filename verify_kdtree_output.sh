#!/bin/bash
set -e

# Send output to log file and terminal
LOG_FILE="verify_kdtree_output.txt"
exec > >(tee -i ${LOG_FILE}) 2>&1

# Ensure we are running from the script directory and starting from a clean state
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
cd "$DIR" || exit 1
rm -rf validation_output
rm -rf ./tests/sage-model-tests/lightcone_3d_SnapNum.png
rm -rf ./tests/sage-model-tests/output_sage_hdf5_one_step_benchmark

./tests/sage-model-tests/run_setup_verify_kdtree_workflow.sh

# 1. Setup Environment
echo "Detecting platform..."
if [[ "$OSTYPE" == "darwin"* ]]; then
    echo "macOS detected - using Homebrew setup"
    source setup_mac.sh
elif [[ -f "/etc/modules" ]] || command -v module >/dev/null 2>&1; then
    echo "HPC/Linux environment detected - assuming modules are loaded or sourcing setup.sh"
    source setup.sh || true
fi

# 2. Define Inputs
TEST_DIR="tests/sage-model-tests"

# Check for benchmark output first (CI default), then standard output (local default)
if [ -d "${TEST_DIR}/output_sage_hdf5_one_step_benchmark/millennium" ]; then
    INPUT_SAGE_DIR="${TEST_DIR}/output_sage_hdf5_one_step_benchmark/millennium"
    echo "Using benchmark SAGE output: ${INPUT_SAGE_DIR}"
elif [ -d "${TEST_DIR}/output_sage_hdf5_one_step/millennium" ]; then
    INPUT_SAGE_DIR="${TEST_DIR}/output_sage_hdf5_one_step/millennium"
    echo "Using standard SAGE output: ${INPUT_SAGE_DIR}"
else
    # Default fallback (will fail check later)
    INPUT_SAGE_DIR="${TEST_DIR}/output_sage_hdf5_one_step_benchmark/millennium"
fi

PARAM_FILE="${TEST_DIR}/input/millennium_sage_hdf5.par"
ALIST_FILE="${TEST_DIR}/input/millennium/trees/millennium.a_list"

# 3. Validation Output Directory
VALIDATION_DIR="validation_output"
rm -rf ${VALIDATION_DIR}
mkdir -p ${VALIDATION_DIR}

echo "Using SAGE input from: ${INPUT_SAGE_DIR}"

if [ ! -f "${INPUT_SAGE_DIR}/model_0.hdf5" ]; then
    echo "ERROR: SAGE output not found at ${INPUT_SAGE_DIR}/model_0.hdf5"
    echo "Please run 'tests/sage-model-tests/run_test_hdf5_one_step_benchmark.sh' (for benchmark output) or 'tests/sage-model-tests/run_test_hdf5_one_step.sh' (for standard output) first to generate input data."
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
    -o ${VALIDATION_DIR}/kdtree_test.h5 \
    --ppc 1000 \
    -v 1
duration_new=$(( SECONDS - start_new ))
echo "NEW sage2kdtree took ${duration_new} seconds."

# Run cli_lightcone to generate lightcone_3d_SnapNum.png for visual verification
echo "========== PHASE 3: Running cli_lightcone =========="
./bin/cli_lightcone \
    --dataset ${VALIDATION_DIR}/kdtree_test.h5 \
    --decmin 0 --decmax 1 --ramin 0 --ramax 1 --zmin 0 --zmax 1 \
    --outdir ${VALIDATION_DIR} --outfile lightcone.h5
    #--outdir ${OUTPUTDIR} --outfields ra dec redshift_cosmological redshift_observed sfr --outfile $RAWNAME-lightcone.h5


# 5: Plot (not benchmarked - visualization only)
echo "========== PHASE 4: Plotting =========="
source .venv/bin/activate
# Field names are case-insensitive in plot_lightcone.py (SnapNum = snapnum)
python3 ./src/plot_lightcone.py ${VALIDATION_DIR}/lightcone.h5 SnapNum
