#!/bin/bash
_F="${BASH_SOURCE[0]:-$0}"
MY_SCRIPTS_DIRECTORY="$(cd "$(dirname "$_F")" && pwd)"
MY_ROOT="$(cd "${MY_SCRIPTS_DIRECTORY}/.." && pwd)"
cd "${MY_ROOT}/tests/sage-model-tests" || exit 1

xctrace record --template 'Time Profiler' --launch -- "${MY_ROOT}/bin/cli_lightcone" \
  --centralgalaxies \
  --dataset output_centralgalaxies_validation/pass1-myhdf5millennium-kdtree.h5 \
  --decmin 0 --decmax 1 --ramin 0 --ramax 1 --zmin 0 --zmax 1 \
  --outdir output_centralgalaxies_validation \
  --outfile profile-lightcone.h5
