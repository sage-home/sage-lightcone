#!/bin/bash

# Platform-aware build script that automatically selects the right setup
# This preserves Linux HPC compatibility while enabling macOS builds

echo "Detecting platform..."

# Detect the platform
if [[ "$OSTYPE" == "darwin"* ]]; then
    echo "macOS detected - using Homebrew setup"
    source setup_mac.sh
    
    # Additional macOS-specific build flags
    export CMAKE_EXTRA_FLAGS="-DCMAKE_OSX_ARCHITECTURES=arm64"
    
elif [[ -f "/etc/modules" ]] || command -v module >/dev/null 2>&1; then
    echo "HPC/Linux environment detected - using module setup"
    source setup.sh
    
    # HPC-specific optimizations - using Debug for debugging
    export CMAKE_EXTRA_FLAGS="-DCMAKE_BUILD_TYPE=Debug"
    
else
    echo "Unknown platform - falling back to basic setup"
    source setup.sh
fi

echo "Platform-specific setup complete"
echo "CMAKE_EXTRA_FLAGS: $CMAKE_EXTRA_FLAGS"

# Clone the selected SAGE repository in place in case it doesn't exist
# That is, clone if doesn't exist, otherwise pull latest

#SAGE_REPO="https://github.com/sage-home/sage-model.git"
SAGE_REPO="https://github.com/sage-home/sage-model-update.git"

# Clone if doesn't exist, otherwise pull latest
if [ ! -d "sage-model" ]; then
    echo "Cloning SAGE from ${SAGE_REPO}..."
    git clone ${SAGE_REPO} sage-model
    if [ $? -ne 0 ]; then
        echo "ERROR: Failed to clone SAGE repository"
        exit 1
    fi
else
    echo "Updating SAGE from ${SAGE_REPO}..."
    cd sage-model
    git pull
    cd ..
fi

# Clean up any old CMake cache files
echo "Cleaning old build artifacts..."
find . -name CMakeCache.txt -exec rm -f {} \;

# Build PugiXML first
echo "Building PugiXML..."
cd dep/pugixml-1.7/scripts
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=.. ${CMAKE_EXTRA_FLAGS} ..
make -j$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)
make install
cd ../../../..

echo "Building main project..."
# CMAKE_BUILD_TYPE = Debug or Release
cmake -B="./bin" -DCMAKE_BUILD_TYPE=Debug \
-DPUGIXML_HOME=./dep/pugixml-1.7 \
${CMAKE_EXTRA_FLAGS} \
.
cd bin
make clean
make -j$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)

# Build SAGE executable from git repository
echo ""
echo "Building SAGE executable..."
cd ..
cd sage-model

# Set environment based on platform detected earlier
if [[ "$OSTYPE" == "darwin"* ]]; then
    echo "Building SAGE for macOS..."

    # Find HDF5 directory (try multiple common locations)
    if command -v brew >/dev/null 2>&1; then
        HDF5_DIR=$(brew --prefix hdf5-mpi 2>/dev/null || brew --prefix hdf5 2>/dev/null || echo "")
    fi

    if [ -z "$HDF5_DIR" ] && [ -d "/opt/homebrew/Cellar/hdf5-mpi" ]; then
        HDF5_DIR="/opt/homebrew"
    elif [ -z "$HDF5_DIR" ] && [ -d "/usr/local/Cellar/hdf5-mpi" ]; then
        HDF5_DIR="/usr/local"
    fi

    if [ -n "$HDF5_DIR" ]; then
        export HDF5_DIR
        echo "Using HDF5_DIR: $HDF5_DIR"
    fi

    # Use system compiler on macOS
    export CC=gcc
else
    echo "Building SAGE for HPC/Linux environment..."
    # Use MPI compiler for HPC environments
    export CC=mpicc
    
    # If HDF5_ROOT is set (common in CI/CMake), map it to HDF5_DIR (used by SAGE Makefile)
    if [ -n "$HDF5_ROOT" ] && [ -z "$HDF5_DIR" ]; then
        export HDF5_DIR="$HDF5_ROOT"
        echo "Mapped HDF5_ROOT to HDF5_DIR: $HDF5_DIR"
    fi
fi

# Build sage
echo "Running make in sage-model directory..."
make

# Check if build was successful
if [ -f "sage" ]; then
    echo "SAGE executable built successfully!"
    echo "SAGE executable location: $(pwd)/sage"

    # Copy to parent bin directory for convenience
    cp sage ../bin/sage
    echo "Copied SAGE executable to bin directory"
else
    echo "ERROR: SAGE executable was not created!"
    exit 1
fi

cd ..

echo ""
echo "Build complete!"
echo "Available executables in bin/:"
ls -la bin/ | grep -E "(cli_lightcone|sage2h5|dstreeinit|sage)$" || echo "No executables found"
