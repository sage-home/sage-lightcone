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

# Build SAGE executable from submodule
echo ""
echo "Building SAGE executable..."
cd ..

# Check if sage-model directory exists
if [ -d "sage-model" ]; then
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
    fi
    
    # Build sage
    echo "Running make all in sage-model directory..."
    make all
    
    # Check if build was successful
    if [ -f "sage" ]; then
        echo "SAGE executable built successfully!"
        echo "SAGE executable location: $(pwd)/sage"
        
        # Copy to parent bin directory for convenience
        cp sage ../bin/sage
        echo "Copied SAGE executable to bin directory"
    else
        echo "Warning: SAGE executable was not created!"
    fi
    
    cd ..
else
    echo "Warning: sage-model directory not found. Skipping SAGE build."
    echo "To build SAGE, initialize the submodule with: git submodule update --init --recursive"
fi

echo ""
echo "Build complete!"
echo "Available executables in bin/:"
ls -la bin/ | grep -E "(cli_lightcone|sage2h5|dstreeinit|sage)$" || echo "No executables found"