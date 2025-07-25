#!/bin/bash

echo "Compiling h5_to_csv with parallel HDF5 support..."

# Source the environment
source ../../setup.sh

# First try with h5pcc (parallel HDF5 compiler wrapper)
if command -v h5pcc &> /dev/null; then
    echo "Trying with h5pcc (parallel HDF5 compiler)..."
    h5pcc -std=c++17 -O2 -Wall -Wextra -o h5_to_csv h5_to_csv.cpp
    if [ $? -eq 0 ]; then
        echo "Compilation with h5pcc successful!"
        echo "You can now run: ./h5_to_csv input.h5 output.csv"
        exit 0
    fi
fi

# Try with mpicxx and explicit HDF5 paths
echo "Trying with mpicxx and explicit HDF5 paths..."
mpicxx -std=c++17 -O2 -Wall -Wextra \
    -I${HDF5_DIR}/include \
    -L${HDF5_DIR}/lib \
    -lhdf5 -lmpi \
    -o h5_to_csv h5_to_csv.cpp

if [ $? -eq 0 ]; then
        cp h5_to_csv ../../runtime/bin
        echo "Compilation with mpicxx successful!"
        echo "h5_to_csv installed to ../../runtime/bin"
else
    echo "Compilation failed. Let's try a different approach..."
    
    # Try with mpicc (C compiler) but as C++
    echo "Trying with mpicc..."
    mpicc -x c++ -std=c++17 -O2 -Wall -Wextra \
        -I${HDF5_DIR}/include \
        -L${HDF5_DIR}/lib \
        -lhdf5 -lstdc++ \
        -o h5_to_csv h5_to_csv.cpp
    
    if [ $? -eq 0 ]; then
        cp h5_to_csv ../../runtime/bin
        echo "Compilation with mpicc successful!"
        echo "h5_to_csv installed to ../../runtime/bin"
    else
        echo "All compilation attempts failed."
        echo "Please check that HDF5 and MPI are properly installed and configured."
        echo "You may need to adjust the HDF5_DIR environment variable."
    fi
fi
