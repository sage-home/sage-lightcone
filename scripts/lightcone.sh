#!/bin/bash
module purge &> /dev/null
module load gcc/12.2.0
module load boost/1.81.0
module load openmpi/4.1.4
module load hdf5/1.14.0
module load gcccore/12.2.0
module load gsl/2.7

export MY_INVOKING_SCRIPT=${BASH_SOURCE:-$0}
export MY_INVOKING_SCRIPTS_DIRECTORY=$(dirname $MY_INVOKING_SCRIPT)
export MY_RUNTIME_DIRECTORY=$(cd ${MY_INVOKING_SCRIPTS_DIRECTORY} && pwd)

${MY_RUNTIME_DIRECTORY}/bin/cli_lightcone $*
