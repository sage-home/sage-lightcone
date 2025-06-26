#!/bin/bash
module purge &> /dev/null
module load gcc/11.3.0
module load boost/1.79.0
module load openmpi/4.1.4
module load soci/4.0.3
module load gsl/2.7
module load hdf5/1.13.1
module load gcccore/12.2.0
module load cfitsio/4.2.0
#module load postgresql/14.4
module load cmake/3.24.3
module load git-lfs/3.2.0
export LD_LIBRARY_PATH=/apps/modules/software/ICU/71.1-GCCcore-11.3.0/lib:$LD_LIBRARY_PATH
