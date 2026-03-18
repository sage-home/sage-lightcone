#!/bin/bash

# What's my code's root directory
export MY_SCRIPT=${BASH_SOURCE:-$0}
export MY_SCRIPTS_DIRECTORY=$(dirname $MY_SCRIPT)
export MY_ROOT=$(cd ${MY_SCRIPTS_DIRECTORY} && pwd)

# Check for Python virtual environment and install if needed
if [[ ! -d "$MY_ROOT/.venv" ]]; then
    echo "Virtual environment not found. Running install.sh..."
    if [[ -f "$MY_ROOT/install.sh" ]]; then
        cd "$MY_ROOT"
        chmod +x install.sh
        ./install.sh
        if [[ $? -eq 0 ]]; then
            echo "✅ Virtual environment created successfully"
        else
            echo "❌ Failed to create virtual environment"
            return 1
        fi
    else
        echo "❌ install.sh not found. Cannot create virtual environment."
        return 1
    fi
else
    echo "✅ Virtual environment found at $MY_ROOT/.venv"
fi

# Activate virtual environment
source $MY_ROOT/.venv/bin/activate
echo "Python virtual environment activated"
pip3 list

# Mac-specific setup script for M3 Mac with Homebrew
# This replaces the module-based setup.sh for macOS development

echo "Setting up environment for M3 Mac..."

# Ensure Homebrew is in PATH
export PATH="/opt/homebrew/bin:/opt/homebrew/sbin:$PATH"

# Set compiler preferences for Apple Silicon
export CC="clang"
export CXX="clang++"

# HDF5-MPI paths - robust discovery
if [[ -z "$USE_MPI" ]]; then
    USE_MPI="no"
fi

# Only look for parallel HDF5 if USE_MPI is "yes"
if [[ "$USE_MPI" == "yes" ]]; then
    HDF5_CANDIDATES=(
        "/opt/homebrew/Cellar/hdf5-mpi"         # Homebrew hdf5-mpi (any version)
        "/opt/homebrew"                         # Homebrew standard location
        "/usr/local"                            # Alternative location
    )

    HDF5_ROOT=""
    for candidate in "${HDF5_CANDIDATES[@]}"; do
        if [[ "$candidate" == "/opt/homebrew/Cellar/hdf5-mpi" ]]; then
            # Find the latest version in the Cellar
            if [[ -d "$candidate" ]]; then
                latest_version=$(ls "$candidate" | sort -V | tail -1)
                if [[ -n "$latest_version" && -d "$candidate/$latest_version" ]]; then
                    HDF5_ROOT="$candidate/$latest_version"
                    echo "Found HDF5-MPI at: $HDF5_ROOT"
                    break
                fi
            fi
        elif [[ -f "$candidate/bin/h5cc" ]]; then
            # Check if this is actually parallel HDF5 (has h5pcc)
            if [[ -f "$candidate/bin/h5pcc" ]]; then
                HDF5_ROOT="$candidate"
                echo "Found HDF5-MPI at: $HDF5_ROOT"
                break
            fi
        fi
    done

    if [[ -n "$HDF5_ROOT" ]]; then
        export HDF5_ROOT
        export HDF5_DIR="$HDF5_ROOT"
        export PATH="$HDF5_ROOT/bin:$PATH"
    else
        echo "Warning: HDF5-MPI not found. Using system PATH."
        if ! command -v h5pcc &> /dev/null; then
             # Only warn if we really wanted MPI
             echo "Warning: h5pcc not found. Please install HDF5 with MPI support: brew install hdf5-mpi"
        fi
    fi
else
    echo "Skipping HDF5-MPI discovery (USE_MPI=no)"
fi

# OpenMPI paths  
export MPI_ROOT="/opt/homebrew"
export PATH="$MPI_ROOT/bin:$PATH"

# Boost paths
export BOOST_ROOT="/opt/homebrew"

# CMake paths - robust discovery
CMAKE_CANDIDATES=(
    "/opt/homebrew"                          # Homebrew on Apple Silicon
    "/usr/local"                            # Homebrew on Intel
    "/Applications/CMake.app/Contents"      # CMake.app installation
    "/usr"                                  # System installation
)

CMAKE_ROOT=""
for candidate in "${CMAKE_CANDIDATES[@]}"; do
    if [[ -x "$candidate/bin/cmake" ]]; then
        CMAKE_ROOT="$candidate"
        echo "Found CMake at: $CMAKE_ROOT"
        break
    fi
done

if [[ -n "$CMAKE_ROOT" ]]; then
    export CMAKE_ROOT
    export PATH="$CMAKE_ROOT/bin:$PATH"
else
    echo "Warning: CMake not found in expected locations. Using system PATH."
    # Check if cmake is available in PATH
    if ! command -v cmake &> /dev/null; then
        echo "Error: cmake not found in PATH either. Please install CMake."
        echo "You can install it via:"
        echo "  - Homebrew: brew install cmake"
        echo "  - Download from: https://cmake.org/download/"
        return 1
    fi
fi

# Set PKG_CONFIG_PATH to help CMake find libraries
export PKG_CONFIG_PATH="/opt/homebrew/lib/pkgconfig:$PKG_CONFIG_PATH"

# Set library and include paths
# IMPORTANT: Only add generic paths if NOT doing a serial build, 
# as they might leak MPI headers from system/brew into the build
if [[ "$USE_MPI" == "yes" ]]; then
    export LIBRARY_PATH="/opt/homebrew/lib:$LIBRARY_PATH"
    export CPATH="/opt/homebrew/include:$CPATH"
fi

echo "Environment configured for:"
if [[ "$USE_MPI" != "yes" ]]; then
    echo "  - MPI: Disabled (USE_MPI=no)"
    # Check for serial HDF5
    H5CC_SERIAL=$(brew --prefix hdf5 2>/dev/null)/bin/h5cc
    H5CC_PARALLEL=$(brew --prefix hdf5-mpi 2>/dev/null)/bin/h5pcc
    
    if [[ -x "$H5CC_SERIAL" ]]; then
         echo "  - HDF5 (Serial): $($H5CC_SERIAL -showconfig 2>/dev/null | grep 'HDF5 Version' || echo 'Found')"
    elif [[ -x "$H5CC_PARALLEL" ]]; then
         echo "  - HDF5 (Parallel-as-Serial): $($H5CC_PARALLEL -showconfig 2>/dev/null | grep 'HDF5 Version' || echo 'Found')"
    else
         echo "  - HDF5 (Serial/Parallel): Not found in standard brew location"
    fi
else
    echo "  - HDF5-MPI: $(h5cc -showconfig 2>/dev/null | grep 'HDF5 Version' || echo 'Not found')"
    echo "  - OpenMPI: $(mpicc --version 2>/dev/null | head -1 || echo 'Not found')"
fi
echo "  - Boost: $(brew list boost --versions 2>/dev/null || echo 'Not found')"
echo "  - CMake: $(cmake --version 2>/dev/null | head -1 || echo 'Not found')"

# Validate critical dependencies
validation_failed=false

if ! command -v cmake &> /dev/null; then
    echo "❌ CMake validation failed"
    validation_failed=true
else
    echo "✅ CMake: $(cmake --version | head -1)"
fi

if [[ "$USE_MPI" == "yes" ]]; then
    if ! command -v h5cc &> /dev/null; then
        echo "❌ HDF5 validation failed"
        validation_failed=true
    else
        # Check if HDF5 has parallel support
        if h5cc -showconfig 2>/dev/null | grep -i "Parallel HDF5.*ON"; then
            echo "✅ HDF5-MPI: Parallel support enabled"
        else
            echo "⚠️  HDF5: Found but parallel support may not be enabled"
        fi
    fi

    if ! command -v mpicc &> /dev/null; then
        echo "❌ MPI validation failed"
        validation_failed=true
    else
        echo "✅ MPI: $(mpicc --version 2>/dev/null | head -1)"
    fi
else
    echo "✅ MPI: Disabled explicitly"
    
    # Check for Serial HDF5
    H5CC_SERIAL=$(brew --prefix hdf5 2>/dev/null)/bin/h5cc
    if [[ -x "$H5CC_SERIAL" ]]; then
        echo "✅ HDF5 (Serial): Found at $(dirname $(dirname $H5CC_SERIAL))"
        # Export HDF5_ROOT for scripts that rely on setup_mac.sh
        export HDF5_ROOT=$(dirname $(dirname $H5CC_SERIAL))
        export HDF5_DIR="$HDF5_ROOT"
    else
        # Fallback for "Parallel HDF5 as Serial"
        H5CC_PARALLEL=$(brew --prefix hdf5-mpi 2>/dev/null)/bin/h5pcc
        if [[ -x "$H5CC_PARALLEL" ]]; then
            echo "✅ HDF5 (Parallel-as-Serial): Found at $(dirname $(dirname $H5CC_PARALLEL))"
            # Export HDF5_ROOT for scripts that rely on setup_mac.sh
            export HDF5_ROOT=$(dirname $(dirname $H5CC_PARALLEL))
            export HDF5_DIR="$HDF5_ROOT"
        else
            echo "❌ HDF5 (Serial) not found! Run 'brew install hdf5' OR 'brew install hdf5-mpi'"
            validation_failed=true
        fi
    fi
fi

if [[ "$validation_failed" == "true" ]]; then
    echo ""
    echo "❌ Some dependencies are missing. Build may fail."
    echo "Consider installing missing dependencies with Homebrew."
    return 1
else
    echo ""
    echo "✅ All critical dependencies found. Ready to build!"
fi
