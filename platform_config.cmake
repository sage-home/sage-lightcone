# Platform-specific CMake configuration helper
# Include this file in your main CMakeLists.txt to get platform-aware settings

# =============================================================================
# macOS Dependency Checks
# =============================================================================
# This section helps users on fresh macOS machines understand what dependencies
# are required and how to install them.
#
# To skip these checks (e.g., if you have deps installed via MacPorts or manually):
#   cmake -DSKIP_MACOS_DEP_CHECK=ON ..

option(SKIP_MACOS_DEP_CHECK "Skip macOS Homebrew dependency checks" OFF)

if(APPLE AND NOT SKIP_MACOS_DEP_CHECK)
    message(STATUS "")
    message(STATUS "Checking macOS dependencies...")

    # Track missing dependencies
    set(MACOS_MISSING_DEPS "")
    set(MACOS_INSTALL_COMMANDS "")

    # --------------------------------------------------------------------------
    # Check 1: Homebrew
    # --------------------------------------------------------------------------
    find_program(HOMEBREW_EXECUTABLE brew)
    if(NOT HOMEBREW_EXECUTABLE)
        message(WARNING
            "Homebrew not found! Homebrew is required to install dependencies on macOS.\n"
            "Install Homebrew with:\n"
            "  /bin/bash -c \"$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)\"\n"
            "Then re-run cmake."
        )
        message(FATAL_ERROR "Cannot continue without Homebrew on macOS.")
    else()
        message(STATUS "  [OK] Homebrew found: ${HOMEBREW_EXECUTABLE}")
    endif()

    # --------------------------------------------------------------------------
    # Helper function to check Homebrew packages
    # --------------------------------------------------------------------------
    function(check_homebrew_package package_name formula_name install_note)
        execute_process(
            COMMAND brew list ${formula_name}
            RESULT_VARIABLE BREW_RESULT
            OUTPUT_QUIET
            ERROR_QUIET
        )
        if(NOT BREW_RESULT EQUAL 0)
            message(STATUS "  [MISSING] ${package_name} (${formula_name})")
            set(MACOS_MISSING_DEPS "${MACOS_MISSING_DEPS};${package_name}" PARENT_SCOPE)
            if(install_note)
                set(MACOS_INSTALL_COMMANDS "${MACOS_INSTALL_COMMANDS}\n  brew install ${formula_name}  # ${install_note}" PARENT_SCOPE)
            else()
                set(MACOS_INSTALL_COMMANDS "${MACOS_INSTALL_COMMANDS}\n  brew install ${formula_name}" PARENT_SCOPE)
            endif()
        else()
            # Get version info
            execute_process(
                COMMAND brew info ${formula_name} --json=v2
                OUTPUT_VARIABLE BREW_INFO
                OUTPUT_STRIP_TRAILING_WHITESPACE
                ERROR_QUIET
            )
            message(STATUS "  [OK] ${package_name}")
        endif()
    endfunction()

    # --------------------------------------------------------------------------
    # Check 2: Required Homebrew packages
    # --------------------------------------------------------------------------

    # HDF5 with MPI support (critical - must be hdf5-mpi, not plain hdf5)
    execute_process(
        COMMAND brew list hdf5-mpi
        RESULT_VARIABLE HDF5_MPI_RESULT
        OUTPUT_QUIET
        ERROR_QUIET
    )
    execute_process(
        COMMAND brew list hdf5
        RESULT_VARIABLE HDF5_PLAIN_RESULT
        OUTPUT_QUIET
        ERROR_QUIET
    )
    if(NOT HDF5_MPI_RESULT EQUAL 0)
        if(HDF5_PLAIN_RESULT EQUAL 0)
            message(STATUS "  [WRONG] HDF5 found but without MPI support!")
            message(STATUS "          You have 'hdf5' but need 'hdf5-mpi' for parallel I/O.")
            list(APPEND MACOS_MISSING_DEPS "hdf5-mpi")
            set(MACOS_INSTALL_COMMANDS "${MACOS_INSTALL_COMMANDS}\n  brew uninstall hdf5  # Remove non-MPI version first")
            set(MACOS_INSTALL_COMMANDS "${MACOS_INSTALL_COMMANDS}\n  brew install hdf5-mpi  # HDF5 with MPI support (required for parallel I/O)")
        else()
            message(STATUS "  [MISSING] HDF5-MPI")
            list(APPEND MACOS_MISSING_DEPS "hdf5-mpi")
            set(MACOS_INSTALL_COMMANDS "${MACOS_INSTALL_COMMANDS}\n  brew install hdf5-mpi  # HDF5 with MPI support (required for parallel I/O)")
        endif()
    else()
        message(STATUS "  [OK] HDF5-MPI (parallel HDF5)")
    endif()

    # OpenMPI
    check_homebrew_package("OpenMPI" "open-mpi" "MPI implementation for parallel processing")

    # Boost
    check_homebrew_package("Boost" "boost" "C++ libraries for program options, filesystem, etc.")

    # GSL (GNU Scientific Library)
    check_homebrew_package("GSL" "gsl" "GNU Scientific Library for numerical computations")

    # LLVM/Clang (optional but recommended for better C++ support)
    execute_process(
        COMMAND brew list llvm
        RESULT_VARIABLE LLVM_RESULT
        OUTPUT_QUIET
        ERROR_QUIET
    )
    if(NOT LLVM_RESULT EQUAL 0)
        message(STATUS "  [OPTIONAL] LLVM/Clang not installed (recommended for better C++ support)")
        message(STATUS "             Install with: brew install llvm")
    else()
        message(STATUS "  [OK] LLVM/Clang")
    endif()

    # --------------------------------------------------------------------------
    # Check 3: Python virtual environment (for sageimport and plotting)
    # --------------------------------------------------------------------------
    if(NOT EXISTS "${CMAKE_SOURCE_DIR}/.venv")
        message(STATUS "  [OPTIONAL] Python virtual environment not found at .venv")
        message(STATUS "             Create with: python3 -m venv .venv && source .venv/bin/activate")
        message(STATUS "             Then install: pip install numpy h5py matplotlib mpi4py")
    else()
        message(STATUS "  [OK] Python virtual environment (.venv)")
    endif()

    # --------------------------------------------------------------------------
    # Summary: Report missing dependencies
    # --------------------------------------------------------------------------
    list(LENGTH MACOS_MISSING_DEPS NUM_MISSING)
    if(NUM_MISSING GREATER 0)
        message(STATUS "")
        message(WARNING
            "Missing ${NUM_MISSING} required macOS dependencies!\n\n"
            "Install missing dependencies with:${MACOS_INSTALL_COMMANDS}\n\n"
            "After installing, re-run cmake to continue.\n"
            "For a complete setup, run: source setup_mac.sh"
        )
        message(FATAL_ERROR "Cannot continue without required dependencies. See instructions above.")
    else()
        message(STATUS "  All required macOS dependencies found!")
    endif()

    message(STATUS "")
elseif(APPLE AND SKIP_MACOS_DEP_CHECK)
    message(STATUS "Skipping macOS dependency checks (SKIP_MACOS_DEP_CHECK=ON)")
endif()

# =============================================================================
# HPC Environment Detection
# =============================================================================

# Function to detect if we're on an HPC environment
function(detect_hpc_environment result_var)
    # Check for common HPC indicators
    if(DEFINED ENV{SLURM_JOB_ID} OR 
       DEFINED ENV{PBS_JOBID} OR 
       DEFINED ENV{LSB_JOBID} OR
       EXISTS "/etc/slurm/slurm.conf" OR
       EXISTS "/opt/pbs" OR
       DEFINED ENV{MODULEPATH})
        set(${result_var} TRUE PARENT_SCOPE)
    else()
        set(${result_var} FALSE PARENT_SCOPE)
    endif()
endfunction()

# Detect HPC environment
detect_hpc_environment(IS_HPC)

# Set C++ standard based on environment
if(IS_HPC)
    # HPC environments might prefer C++20 for newer features
    set(PLATFORM_CXX_STANDARD 20)
    message(STATUS "HPC environment detected - using C++${PLATFORM_CXX_STANDARD}")
elseif(APPLE)
    # macOS prefers C++17 for better compatibility
    set(PLATFORM_CXX_STANDARD 17)
    message(STATUS "macOS environment detected - using C++${PLATFORM_CXX_STANDARD}")
else()
    # Default to C++17 for general Linux
    set(PLATFORM_CXX_STANDARD 17)
    message(STATUS "Default environment - using C++${PLATFORM_CXX_STANDARD}")
endif()

# Override default if not already set
if(NOT DEFINED CXX_STANDARD_DEFAULT)
    set(CXX_STANDARD_DEFAULT ${PLATFORM_CXX_STANDARD})
endif()