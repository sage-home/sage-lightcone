#!/bin/bash
set -e

# Ensure we are running from the script directory
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
cd "$DIR" || exit 1

rm -rf ./tests/sage-model-tests/lightcone_3d_SnapNum.png
rm -rf ./tests/sage-model-tests/output_sage_hdf5_one_step_benchmark

./tests/sage-model-tests/robust_run_test_hdf5_one_step_benchmark.sh
./verify_kdtree_output.sh
