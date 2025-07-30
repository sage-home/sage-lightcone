#!/bin/bash
if [ "$SYS_ARCH" = "milan" ];then
module purge &> /dev/null
module load gcc/12.2.0
module load boost/1.81.0
module load openmpi/4.1.4
module load hdf5/1.14.0
module load gcccore/12.2.0
module load gsl/2.7
module load cmake/3.24.3
fi
