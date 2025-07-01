#!/bin/bash
source ./setup.sh
# For now the validation tests still exit with 0 
bin/cli_lightcone --dataset tests/dummy.h5 --decmin 100 --decmax 0
bin/cli_lightcone --dataset tests/dummy.h5 --decmin 10 --decmax 0
bin/cli_lightcone --dataset tests/dummy.h5 --decmin -10 --decmax 100

bin/cli_lightcone --dataset tests/dummy.h5 --ramin -10 --ramax 100
bin/cli_lightcone --dataset tests/dummy.h5 --ramin 10 --ramax 0
bin/cli_lightcone --dataset tests/dummy.h5 --ramin -10 --ramax 400
bin/cli_lightcone --dataset tests/dummy.h5 --ramin 0 --ramax 10 --zmin -1 --zmax 200
bin/cli_lightcone --dataset tests/dummy.h5 --ramin 0 --ramax 10 --zmin 3 --zmax 2

bin/cli_lightcone --dataset tests/dummy.h5 --decmin 0 --decmax 10 --ramin 0 --ramax 10 --outfields dummy
