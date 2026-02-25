#!/bin/bash
source ./setup.sh
# For now the validation tests still exit with 0 (cli_lightcone doesn't return non-zero on validation errors yet)
# These tests verify that invalid arguments are detected and reported.
bin/cli_lightcone --dataset tests/dummy.h5 --decmin 100 --decmax 0 > /dev/null
bin/cli_lightcone --dataset tests/dummy.h5 --decmin 10 --decmax 0 > /dev/null
bin/cli_lightcone --dataset tests/dummy.h5 --decmin -10 --decmax 100 > /dev/null

bin/cli_lightcone --dataset tests/dummy.h5 --ramin -10 --ramax 100 > /dev/null
bin/cli_lightcone --dataset tests/dummy.h5 --ramin 10 --ramax 0 > /dev/null
bin/cli_lightcone --dataset tests/dummy.h5 --ramin -10 --ramax 400 > /dev/null
bin/cli_lightcone --dataset tests/dummy.h5 --ramin 0 --ramax 10 --zmin -1 --zmax 200 > /dev/null
bin/cli_lightcone --dataset tests/dummy.h5 --ramin 0 --ramax 10 --zmin 3 --zmax 2 > /dev/null
