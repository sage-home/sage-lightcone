#!/bin/bash

# Platform-aware build script that automatically selects the right setup
# This preserves Linux HPC compatibility while enabling macOS builds

echo "Detecting platform..."

# Detect the platform
if [[ "$OSTYPE" == "darwin"* ]]; then
    echo "macOS detected - using Homebrew setup"
    source setup_mac.sh
    
    # Additional macOS-specific build flags
    CMAKE_EXTRA_FLAGS="-DCMAKE_OSX_ARCHITECTURES=arm64"
    
elif [[ -f "/etc/modules" ]] || command -v module >/dev/null 2>&1; then
    echo "HPC/Linux environment detected - using module setup"
    source setup.sh
    
    # HPC-specific optimizations - using Debug for debugging
    CMAKE_EXTRA_FLAGS="-DCMAKE_BUILD_TYPE=Debug"
    
else
    echo "Unknown platform - falling back to basic setup"
    source setup.sh
fi

# Handle MPI configuration
if [[ -z "$USE_MPI" ]]; then
    USE_MPI="no"
fi

echo "MPI Configuration: USE_MPI=$USE_MPI"

if [[ "$USE_MPI" == "no" ]]; then
    echo "Disabling MPI support..."
    
    # On macOS, handle HDF5 selection (serial vs parallel used as serial)
    if [[ "$OSTYPE" == "darwin"* ]]; then
        SERIAL_HDF5=$(brew --prefix hdf5 2>/dev/null)
        PARALLEL_HDF5=$(brew --prefix hdf5-mpi 2>/dev/null)
        
        if [ -d "$SERIAL_HDF5" ]; then
             echo "Forcing Serial HDF5 from: $SERIAL_HDF5"
             
             # Export standard CMake hints
             export HDF5_ROOT="$SERIAL_HDF5"
             export HDF5_DIR="$SERIAL_HDF5"
             
             # Prepend to PATH/CMAKE_PREFIX_PATH to shadow system HDF5
             export CMAKE_PREFIX_PATH="$SERIAL_HDF5:$CMAKE_PREFIX_PATH"
             
             # Force include/library paths in compiler flags to override system paths
             # This is critical because /opt/homebrew/include (hdf5-mpi) is in default search paths
             export CXXFLAGS="-I$SERIAL_HDF5/include -L$SERIAL_HDF5/lib $CXXFLAGS"
             export CFLAGS="-I$SERIAL_HDF5/include -L$SERIAL_HDF5/lib $CFLAGS"
             
             # Also pass these explicitly to CMake
             CMAKE_HDF5_FLAGS="-DHDF5_ROOT=$SERIAL_HDF5 -DHDF5_NO_SYSTEM_PATHS=ON -DHDF5_USE_STATIC_LIBRARIES=OFF"
        elif [ -d "$PARALLEL_HDF5" ]; then
             echo "Using Parallel HDF5 as Serial from: $PARALLEL_HDF5"
             # If only parallel HDF5 is available, use it but disable MPI in code
             export HDF5_ROOT="$PARALLEL_HDF5"
             CMAKE_HDF5_FLAGS="-DHDF5_ROOT=$PARALLEL_HDF5 -DHDF5_PREFER_PARALLEL=ON"
        fi
    fi

    export CMAKE_EXTRA_FLAGS="${CMAKE_EXTRA_FLAGS} ${CMAKE_HDF5_FLAGS} -DUSE_MPI=OFF -DCMAKE_CXX_FLAGS=\"-DUSE_MPI=0\""
    export USE_MPI_SAGE=""  # Empty string disables MPI in SAGE Makefile (ifdef check)
    
    # Check if we are using parallel HDF5 (which necessitates mpi.h visibility)
    # This applies to both macOS (if using hdf5-mpi) and Linux
    USING_PARALLEL_HDF5=0
    # Check all parallel HDF5 indicators: file paths under HDF5_ROOT, h5pcc in PATH,
    # or an explicit CMake flag. Any one of these is sufficient.
    if [ -f "$HDF5_ROOT/bin/h5pcc" ] || \
       [ -f "$HDF5_ROOT/include/mpi.h" ] || \
       [ -d "$HDF5_ROOT/include/openmpi" ] || \
       command -v h5pcc >/dev/null 2>&1 || \
       [[ "$CMAKE_EXTRA_FLAGS" == *"HDF5_PREFER_PARALLEL=ON"* ]]; then
        USING_PARALLEL_HDF5=1
    fi

    if [ "$USING_PARALLEL_HDF5" -eq 1 ]; then
         echo "Keeping MPI compilers for SAGE (parallel HDF5 headers detected)"
         # Tell cmake to add -DUSE_PARALLEL_HDF5_IN_SERIAL to compile flags.
         # This makes mpi_stub.hh pull in the real mpi.h before defining stubs,
         # avoiding typedef redefinition conflicts when HDF5 later includes mpi.h.
         # Passed as a cmake boolean option (no spaces) to avoid quoting issues.
         export CMAKE_EXTRA_FLAGS="${CMAKE_EXTRA_FLAGS} -DUSE_PARALLEL_HDF5_IN_SERIAL=ON"
         # For SAGE, we must use mpicc if headers include mpi.h
         if command -v mpicc >/dev/null; then
            export CC="mpicc"
            export CXX="mpicxx"
            # Pass compiler via -D flags so cmake uses them even when CMakeLists.txt
            # sets compiler defaults via set() before project() (which would otherwise
            # shadow the CC/CXX env vars).
            export CMAKE_EXTRA_FLAGS="${CMAKE_EXTRA_FLAGS} -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx"
         else
            echo "WARNING: Parallel HDF5 detected but mpicc not found. SAGE build may fail."
            export CC="gcc"
            export CXX="g++"
         fi
    else
        # Serial HDF5 or no HDF5 detected
        if [[ "$OSTYPE" != "darwin"* ]]; then
            export CC="gcc"
            export CXX="g++"
        fi
        # On macOS, we let default compilers stand (clang usually, or gcc from brew)
    fi
else
    echo "Enabling MPI support..."
    export CMAKE_EXTRA_FLAGS="${CMAKE_EXTRA_FLAGS} -DUSE_MPI=ON -DCMAKE_CXX_FLAGS=\"-DUSE_MPI=1\""
    export USE_MPI_SAGE="yes"
    # Ensure MPI compilers are used and pass via -D flags so cmake cache picks them up
    # (a plain set() in CMakeLists.txt before project() would shadow CC/CXX env vars)
    if command -v mpicc >/dev/null; then
        export CC="mpicc"
        export CXX="mpicxx"
        export CMAKE_EXTRA_FLAGS="${CMAKE_EXTRA_FLAGS} -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx"
    fi
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

# Set HDF5 variables for CMake to ensure correct version (serial vs MPI) is picked up
if [[ "$OSTYPE" == "darwin"* ]] && command -v brew >/dev/null 2>&1; then
    if [[ "$USE_MPI" == "no" ]]; then
        HDF5_PREFIX=$(brew --prefix hdf5 2>/dev/null)

        if [ -n "$HDF5_PREFIX" ] && [ -d "$HDF5_PREFIX" ]; then
            echo "Forcing serial HDF5 for CMake: $HDF5_PREFIX"
            export HDF5_ROOT="$HDF5_PREFIX"
            # Explicitly point to serial HDF5 to avoid picking up MPI version
            CMAKE_EXTRA_FLAGS="${CMAKE_EXTRA_FLAGS} -DHDF5_ROOT=$HDF5_PREFIX -DHDF5_PREFER_PARALLEL=OFF"
        fi
        # Parallel-only HDF5 case: already handled above (HDF5_ROOT, mpicc,
        # PREFER_PARALLEL=ON, and -DUSE_MPI=OFF all set in the USE_MPI=no block)
    else
        HDF5_PREFIX=$(brew --prefix hdf5-mpi 2>/dev/null)
        if [ -n "$HDF5_PREFIX" ]; then
            echo "Preferring parallel HDF5 for CMake: $HDF5_PREFIX"
            # Export HDF5_ROOT so FindHDF5 sees it
            export HDF5_ROOT="$HDF5_PREFIX" 
            CMAKE_EXTRA_FLAGS="${CMAKE_EXTRA_FLAGS} -DHDF5_ROOT=$HDF5_PREFIX -DHDF5_PREFER_PARALLEL=TRUE"
        fi
    fi
fi

echo "Building main project..."
# CMAKE_BUILD_TYPE = Debug or Release (defaults to Debug; override via env var)
cmake -B="./bin" -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE:-Debug} \
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

    # Find HDF5 directory for SAGE Makefile (HDF5_DIR variable)
    # Prefer HDF5_ROOT already determined correctly above (avoids brew --prefix
    # returning hypothetical paths for uninstalled packages on newer Homebrew)
    if [ -n "$HDF5_ROOT" ] && [ -d "$HDF5_ROOT" ]; then
        HDF5_DIR="$HDF5_ROOT"
    elif command -v brew >/dev/null 2>&1; then
        if [[ "$USE_MPI" == "no" ]]; then
            HDF5_DIR=$(brew --prefix hdf5 2>/dev/null || brew --prefix hdf5-mpi 2>/dev/null)
        else
            HDF5_DIR=$(brew --prefix hdf5-mpi 2>/dev/null || brew --prefix hdf5 2>/dev/null)
        fi
        HDF5_DIR=${HDF5_DIR:-""}
    fi

    if [ -n "$HDF5_DIR" ]; then
        export HDF5_DIR
        echo "Using HDF5_DIR: $HDF5_DIR"
    fi

    # Use system compiler on macOS unless Parallel HDF5 forces us to use mpicc
    if [ "$USING_PARALLEL_HDF5" -eq 1 ]; then
        echo "SAGE: Parallel HDF5 detected on macOS, forcing mpicc..."
        if command -v mpicc >/dev/null; then
            export CC=mpicc
            export CXX=mpicxx
        else
            echo "WARNING: mpicc not found, falling back to gcc (SAGE build may fail)"
            export CC=gcc
        fi
    else
        export CC=gcc
    fi
else
    echo "Building SAGE for HPC/Linux environment..."
    # Use MPI compiler for HPC environments if enabled
    if [[ "$USE_MPI" == "yes" ]] || [ "$USING_PARALLEL_HDF5" -eq 1 ]; then
       echo "Using MPI compiler (enabled=$USE_MPI, parallel_hdf5=$USING_PARALLEL_HDF5)..."
       export CC=mpicc
    else
       export CC=gcc
    fi
    
    # If HDF5_ROOT is set (common in CI/CMake), map it to HDF5_DIR (used by SAGE Makefile)
    if [ -n "$HDF5_ROOT" ] && [ -z "$HDF5_DIR" ]; then
        export HDF5_DIR="$HDF5_ROOT"
        echo "Mapped HDF5_ROOT to HDF5_DIR: $HDF5_DIR"
    fi
fi

# Build sage
# Pass CC explicitly on the make command line so it overrides any CC=gcc that the
# SAGE Makefile sets internally when USE-MPI is empty (env-exported CC is not enough
# because Makefile variable assignments take precedence over the environment).
echo "Running make in sage-model directory..."
echo "Using CC='$CC' USE-MPI='$USE_MPI_SAGE'"
make CC="$CC" USE-MPI="$USE_MPI_SAGE"

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
