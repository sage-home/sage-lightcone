# Platform-specific CMake configuration helper
# Include this file in your main CMakeLists.txt to get platform-aware settings

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