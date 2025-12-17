#!/bin/bash

# Validation script for sage2kdtree
# Compares sage2kdtree output against:
# 1. HDF5 workflow baseline (sageh5toh5 → sageimport → dstreeinit)
# 2. Binary workflow baseline (sage2h5 → sageimport → dstreeinit) - AUTHORITATIVE

export MY_SCRIPT=${BASH_SOURCE:-$0}
export MY_SCRIPTS_DIRECTORY=$(dirname $MY_SCRIPT)
export MY_ROOT=$(cd ${MY_SCRIPTS_DIRECTORY}/../.. && pwd)
export MY_SCRIPT=
export MY_SCRIPTS_DIRECTORY=

if [[ "$OSTYPE" == "darwin"* ]]; then
    echo "macOS detected - using Homebrew setup"
    source ${MY_ROOT}/setup_mac.sh
elif [[ -f "/etc/modules" ]] || command -v module >/dev/null 2>&1; then
    echo "HPC/Linux environment detected - using module setup"
    source ${MY_ROOT}/setup.sh
else
    echo "Unknown platform - falling back to basic setup"
    source ${MY_ROOT}/setup.sh
fi

echo "==================================================================="
echo "Validating sage2kdtree against baselines"
echo "==================================================================="

SAGE2KDTREE_OUTPUT="output_sage_hdf5/myhdf5millennium-kdtree-onestep.h5"

# Check sage2kdtree output exists
if [ ! -f "$SAGE2KDTREE_OUTPUT" ]; then
    echo "ERROR: sage2kdtree output not found: $SAGE2KDTREE_OUTPUT"
    echo "Run ./run_test_hdf5_one_step.sh first"
    exit 1
fi

echo ""
echo "==================================================================="
echo "BASELINE 1: HDF5 Workflow (sageh5toh5 → sageimport → dstreeinit)"
echo "==================================================================="

# Run HDF5 workflow baseline
echo "Running HDF5 baseline workflow..."
export RAWNAME=myhdf5millennium
export OUTPUTDIR=output_sage_hdf5_baseline

rm -rf ${OUTPUTDIR}
mkdir -p ${OUTPUTDIR}

# Copy SAGE HDF5 output from main test
cp -r output_sage_hdf5/millennium ${OUTPUTDIR}/

# Run old HDF5 workflow
${MY_ROOT}/bin/sageh5toh5 -m convert -v 2 -s ${OUTPUTDIR}/millennium -p input/millennium_sage_hdf5.par -a input/millennium/trees/millennium.a_list -o ${OUTPUTDIR}/${RAWNAME}-depthfirstordered.h5

${MY_ROOT}/bin/sageimport --settings ${OUTPUTDIR}/${RAWNAME}-depthfirstordered_import_settings.xml

${MY_ROOT}/bin/dstreeinit --mode kdtree --tree ${OUTPUTDIR}/out_${RAWNAME}-depthfirstordered.h5 --sage ${OUTPUTDIR}/${RAWNAME}-bysnap.h5 --output ${OUTPUTDIR}/${RAWNAME}-kdtree.h5

HDF5_BASELINE="${OUTPUTDIR}/${RAWNAME}-kdtree.h5"

if [ ! -f "$HDF5_BASELINE" ]; then
    echo "ERROR: HDF5 baseline not created"
    exit 1
fi

echo ""
echo "Comparing sage2kdtree vs HDF5 baseline:"
echo "  sage2kdtree: $SAGE2KDTREE_OUTPUT"
echo "  HDF5 base:   $HDF5_BASELINE"
echo ""

h5diff -v "$SAGE2KDTREE_OUTPUT" "$HDF5_BASELINE" 2>&1 | head -100

if [ ${PIPESTATUS[0]} -eq 0 ]; then
    echo ""
    echo "✅ PASS: sage2kdtree matches HDF5 workflow baseline!"
else
    echo ""
    echo "❌ FAIL: sage2kdtree differs from HDF5 workflow baseline"
    echo "See differences above"
fi

echo ""
echo "==================================================================="
echo "BASELINE 2: Binary Workflow (AUTHORITATIVE)"
echo "==================================================================="
echo "Running binary baseline workflow (this is the authoritative test)..."

# Run binary workflow baseline
./run_test_binary.sh > /dev/null 2>&1

BINARY_BASELINE="output_sage_binary/mybinarymillennium-kdtree.h5"

if [ ! -f "$BINARY_BASELINE" ]; then
    echo "ERROR: Binary baseline not created"
    exit 1
fi

echo ""
echo "Comparing sage2kdtree vs Binary baseline (AUTHORITATIVE):"
echo "  sage2kdtree:  $SAGE2KDTREE_OUTPUT"
echo "  Binary base:  $BINARY_BASELINE"
echo ""

h5diff -v "$SAGE2KDTREE_OUTPUT" "$BINARY_BASELINE" 2>&1 | head -100

if [ ${PIPESTATUS[0]} -eq 0 ]; then
    echo ""
    echo "✅✅ PASS: sage2kdtree matches AUTHORITATIVE binary workflow baseline!"
    echo ""
    echo "==================================================================="
    echo "SUCCESS: All validations passed!"
    echo "==================================================================="
else
    echo ""
    echo "❌ FAIL: sage2kdtree differs from AUTHORITATIVE binary workflow baseline"
    echo "See differences above"
    exit 1
fi
