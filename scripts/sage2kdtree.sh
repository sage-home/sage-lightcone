#!/bin/bash

export MY_INVOKING_SCRIPT=${BASH_SOURCE:-$0}
export MY_INVOKING_SCRIPTS_DIRECTORY=$(cd "$(dirname "$MY_INVOKING_SCRIPT")" && pwd)
export MY_ROOT=$(cd "${MY_INVOKING_SCRIPTS_DIRECTORY}/.." && pwd)

if [[ "$OSTYPE" == "darwin"* ]]; then
    source "${MY_ROOT}/setup_mac.sh"
else
    module purge &> /dev/null
    module load gcc/12.2.0
    module load boost/1.81.0
    module load openmpi/4.1.4
    module load hdf5/1.14.0
    module load gcccore/12.2.0
    module load gsl/2.7   # required by SAGE (sage-model submodule)
fi

"${MY_ROOT}/bin/sage2kdtree" "$@"
