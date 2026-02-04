#!/bin/bash
set -e

# Send output to log file and terminal
LOG_FILE="verify_kdtree_output.txt"
exec > >(tee -i ${LOG_FILE}) 2>&1

# Validation script for sage2kdtree refactoring
# Compares output of new 'sage2kdtree' vs legacy-copy 'sage2kdtree_four'

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
INPUT_SAGE_DIR="${TEST_DIR}/output_sage_hdf5_one_step_benchmark/millennium"
PARAM_FILE="${TEST_DIR}/input/millennium_sage_hdf5.par"
ALIST_FILE="${TEST_DIR}/input/millennium/trees/millennium.a_list"

# 3. Validation Output Directory
VALIDATION_DIR="validation_output"
rm -rf ${VALIDATION_DIR}
mkdir -p ${VALIDATION_DIR}

echo "Using SAGE input from: ${INPUT_SAGE_DIR}"

if [ ! -f "${INPUT_SAGE_DIR}/model_0.hdf5" ]; then
    echo "ERROR: SAGE output not found at ${INPUT_SAGE_DIR}/model_0.hdf5"
    echo "Please run 'tests/sage-model-tests/run_test_hdf5_one_step.sh' first to generate input data."
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

# 5. Run Legacy Code (sage2kdtree_four)
echo "=========================================================="
echo "Running LEGACY sage2kdtree_four (the baseline)..."
echo "=========================================================="
start_legacy=$SECONDS
./bin/sage2kdtree_four \
    -s ${INPUT_SAGE_DIR} \
    -p ${PARAM_FILE} \
    -a ${ALIST_FILE} \
    -o ${VALIDATION_DIR}/kdtree_baseline.h5 \
    --ppc 1000 \
    -v 1
duration_legacy=$(( SECONDS - start_legacy ))
echo "LEGACY sage2kdtree_four took ${duration_legacy} seconds."

# 6. Compare Outputs
echo ""
echo "=========================================================="
echo "Comparing outputs with python script..."
echo "=========================================================="
echo "Timings Summary:"
echo "  NEW (sage2kdtree):       ${duration_new}s"
echo "  LEGACY (sage2kdtree_four): ${duration_legacy}s"
echo "=========================================================="

source .venv/bin/activate
python3 verify_kdtree.py ${VALIDATION_DIR}/kdtree_baseline.h5 ${VALIDATION_DIR}/kdtree_test.h5
