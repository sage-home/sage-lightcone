#!/bin/bash
#bin/cli_lightcone --help
#bin/cli_lightcone --notknow
bin/cli_lightcone --dataset tests/dummy.h5 --decmin 100 --decmax 0
bin/cli_lightcone --dataset tests/dummy.h5 --decmin 10 --decmax 0
bin/cli_lightcone --dataset tests/dummy.h5 --decmin -10 --decmax 100

bin/cli_lightcone --dataset tests/dummy.h5 --ramin -10 --ramax 100
bin/cli_lightcone --dataset tests/dummy.h5 --ramin 10 --ramax 0
bin/cli_lightcone --dataset tests/dummy.h5 --ramin -10 --ramax 400
bin/cli_lightcone --dataset tests/dummy.h5 --ramin 0 --ramax 10 --zmin -1 --zmax 200
bin/cli_lightcone --dataset tests/dummy.h5 --ramin 0 --ramax 10 --zmin 3 --zmax 2

bin/cli_lightcone --dataset tests/dummy.h5 --decmin 0 --decmax 10 --ramin 0 --ramax 10 --outfields dummy
#bin/cli_lightcone --dataset /fred/oz114/kdtao/bigdata/millennium-public-sage-ozstar-sed-kdtree.h5 \
#--outdir 16x1 \
#--outfile 16x1 \
#--unique false --decmin 0 --decmax 10 --ramin 0 --ramax 10 --zmin 0 --zmax 1
