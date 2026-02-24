# This is the CMakeList.txt file used to build the project
#    - We use this old version because it is the  -
#    -   latest version availale on Travis CI     -
cmake_minimum_required(VERSION 3.10)

# Set CMake policies to suppress warnings
# CMP0074: find_package uses <PackageName>_ROOT variables (CMake 3.12+)
if(POLICY CMP0074)
    cmake_policy(SET CMP0074 NEW)
endif()
# CMP0144: find_package uses upper-case <PACKAGENAME>_ROOT variables (CMake 3.27+)
if(POLICY CMP0144)
    cmake_policy(SET CMP0144 NEW)
endif()
# CMP0167: FindBoost module is deprecated, use BoostConfig.cmake (CMake 3.30+)
if(POLICY CMP0167)
    cmake_policy(SET CMP0167 NEW)
endif()

# Read the file with the macro definitions which specify
# the behaviour specific to this project and then run
# the initialization macro
include(${CMAKE_SOURCE_DIR}/project.cmake)

# Include platform-specific configuration
include(${CMAKE_SOURCE_DIR}/platform_config.cmake)

# Add the 'cmake' directory in the project directory to the module path
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_SOURCE_DIR};${CMAKE_SOURCE_DIR}/.cmake")

# Load and run a macro to validate the initialization
include ("init_validate")
init_validate()

# Set install directory to ./runtime
# (user can still override with DESTDIR with make to set installation directory)
set( CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/runtime" )

# Set output directory for executables to the build directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# Set default compilers (must come before project line)
# Use CACHE without FORCE so that -DCMAKE_C_COMPILER=mpicc on the command line
# (or in CMAKE_EXTRA_FLAGS) takes precedence over these defaults.
# A plain set() before project() would shadow the cache and override -D flags.
set(CMAKE_C_COMPILER   "gcc" CACHE STRING "C compiler")
set(CMAKE_CXX_COMPILER "g++" CACHE STRING "CXX compiler")

# For Clang build
#set(CMAKE_C_COMPILER   "clang" CACHE STRING "C compiler")
#set(CMAKE_CXX_COMPILER "clang++" CACHE STRING "CXX compiler")

# Compiler flags - separate C and C++ flags
# Note: -march=native and -mtune=native can cause issues on Apple Silicon
if(CMAKE_SYSTEM_PROCESSOR MATCHES "arm64|aarch64")
    # Apple Silicon specific flags
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu99 -Wextra")
else()
    # Intel/AMD specific flags
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu99 -Wextra -mtune=native -march=native")
endif()
#set(CMAKE_CXX_FLAGS "-std=c++14 -Wextra -mtune=native -march=native")
#set(CMAKE_CXX_FLAGS "-std=c++17 -Wextra -mtune=native -march=native")
#set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_EXTENSIONS OFF)
# Configure the cmake C++ compile features
# To change the C++ standard at build time (to c++23, for example), run cmake like so:
#    `cmake -DCMAKE_CXX_STANDARD=23 .."
set(CXX_STANDARD_DEFAULT 17)
set(CMAKE_CXX_STANDARD ${CXX_STANDARD_DEFAULT} CACHE STRING "The default C++ standard used")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if (NOT CMAKE_CXX_STANDARD EQUAL CXX_STANDARD_DEFAULT)
    message(STATUS "Selected C++ standard: c++" ${CMAKE_CXX_STANDARD})
else()
    message(STATUS "Using default C++ standard: c++" ${CXX_STANDARD_DEFAULT})
endif()
# Set the project name
set(GBP_PROJECT_NAME "cli_lightcone")

# Declare the project name
project (${GBP_PROJECT_NAME})

# Check for minimum versions of compilers to avoid ABI conflicts
message(STATUS "Checking compiler:")
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(VERSION_MIN "5.0.0")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS ${VERSION_MIN})
        message(FATAL_ERROR "Insufficient gcc version (${CMAKE_CXX_COMPILER_VERSION}<${VERSION_MIN})")
    endif()
elseif((CMAKE_CXX_COMPILER_ID STREQUAL "Clang") OR (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang"))
    set(VERSION_MIN "3.9.0")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS ${VERSION_MIN})
        message(FATAL_ERROR "Insufficient Clang version (${CMAKE_CXX_COMPILER_VERSION}<${VERSION_MIN})")
    endif()
else()
    message(FATAL_ERROR "Unsupported compiler or platform (${CMAKE_CXX_COMPILER_ID})")
endif()
message(STATUS "   ID     : ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "   VERSION: ${CMAKE_CXX_COMPILER_VERSION}")

# Set the backend to kdtrees
message (STATUS "Setting backend to use KDTrees")
add_definitions(-DBACKEND_KDTREE)

# Set no progress for now until progress works correctly
# Currently progress intermittently waits forever for some threads to join due to broken MPI passing being broken
# TODO: Fixme
add_definitions("-DNO_PROGRESS")
if(USE_MPI)
    add_definitions("-DPARALLELHDF5")
    add_definitions("-DUSE_MPI")
endif()
# When building without MPI but against parallel HDF5 headers (which include mpi.h),
# mpi_stub.hh needs this define so it pulls in the real mpi.h before defining stubs,
# preventing typedef redefinition conflicts.
if(USE_PARALLEL_HDF5_IN_SERIAL)
    add_definitions("-DUSE_PARALLEL_HDF5_IN_SERIAL")
endif()

# Set no checkpoint for now as it is broken in testing
# TODO: Fixme
add_definitions("-DNO_CHECKPOINT")

# Set the debug/release defines
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Compiling in DEBUG mode...")
    # RS dont want big log files actually
    #add_definitions(-DNINSTRUMENT -DNSTACKTRACE -DNMEMDEBUG -DNMEMOPS -DNMEMSTATS)
    add_definitions(-DNLOGDEBUG -DNLOGTRIVIAL -DNINSTRUMENT -DNSTACKTRACE -DNMEMDEBUG -DNMEMOPS -DNMEMSTATS)
    # Debug compiler flags
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -g")
    # Suppress duplicate library warnings on macOS
    if(APPLE)
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-no_warn_duplicate_libraries")
    endif()
else()
    message(STATUS "Compiling in RELEASE mode...")
    add_definitions(-DNDEBUG -DNLOGDEBUG -DNLOGTRIVIAL -DNINSTRUMENT -DNSTACKTRACE -DNMEMDEBUG -DNMEMOPS -DNMEMSTATS)

    # CURRENTLY NDEBUG BREAKS THE PROJECT AT RUNTIME
    # TODO: Remove this once we find out why NDEBUG breaks at runtime
    #message(FATAL_ERROR "COMPILING IN RELEASE MODE IS KNOWN TO BREAK AT RUNTIME. -DNDEBUG IS BROKEN.")

    # Release compiler flags
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O1")
    # Suppress duplicate library warnings on macOS
    if(APPLE)
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-no_warn_duplicate_libraries")
    endif()
endif()

# Load all macros and functions
message(STATUS "Loading macros & functions...")
include ("build")          # Main build macros
include ("3rd_party")      # 3rd party library support functions
include ("custom_targets") # Custom-target macros
message(STATUS "Done.")


# Generate compile commands.json file
# (needed by things like clang - tidy)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ===========================================
# Execute all project specific configurations
# (set options, library dependencies, etc.)
# ===========================================
project_options_and_dependencies()
# ===========================================

# Initialize &/or build any needed external libraries.
process_externs( ${CMAKE_SOURCE_DIR} )

# ============= Process the project =============
# (look in cmake/build.cmake for these functions)

# Generate the list of directories where header files
# are located, as well as a list of all header files
process_headers( ${CMAKE_SOURCE_DIR} )

# Generate all library and executable targets
process_targets( ${CMAKE_SOURCE_DIR} )

# ============= Add custom targets ==============
# (look in cmake/custom_targets.cmake for these functions)

# Print status message
message(STATUS "")
message(STATUS "Adding custom targets...")

# Add clang-format, clang-tidy, etc
add_custom_clang()

# Add documentation build
add_custom_docs()

# ============= Add SAGE executable target =============
# Build the SAGE executable from the sage-model submodule
message(STATUS "Adding SAGE executable target...")

# Check if sage-model directory exists
if(EXISTS "${CMAKE_SOURCE_DIR}/sage-model")
    # Check if sage-model submodule is properly initialized
    if(NOT EXISTS "${CMAKE_SOURCE_DIR}/sage-model/Makefile")
        message(WARNING "sage-model directory exists but appears to be empty or incomplete!")
        message(WARNING "This usually means the Git submodule wasn't properly initialized.")
        message(WARNING "Please run: git submodule update --init --recursive")
        message(WARNING "SAGE executable will not be built.")
    else()
        # Determine platform-specific make command and environment
    if(${PLATFORM_MACOS})
        # On macOS, we need to ensure HDF5 is found properly
        set(SAGE_ENV_VARS 
            "CC=${CMAKE_C_COMPILER}"
            "HDF5_DIR=${HDF5_DIR}"
        )
        if(CMAKE_BUILD_TYPE STREQUAL "Debug")
            list(APPEND SAGE_ENV_VARS "CCFLAGS=-g -O0")
        else()
            list(APPEND SAGE_ENV_VARS "CCFLAGS=-O2")
        endif()
    else()
        # For HPC environments, use mpicc if available
        set(SAGE_ENV_VARS 
            "CC=mpicc"
        )
        if(CMAKE_BUILD_TYPE STREQUAL "Debug")
            list(APPEND SAGE_ENV_VARS "CCFLAGS=-g -O0")
        else()
            list(APPEND SAGE_ENV_VARS "CCFLAGS=-O2 -march=native")
        endif()
    endif()

    # Create the sage executable target
    add_custom_target(sage_executable
        COMMAND ${CMAKE_COMMAND} -E env ${SAGE_ENV_VARS} make all
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/sage-model
        COMMENT "Building SAGE executable from sage-model submodule"
        VERBATIM
    )

    # Make sage_executable part of the default build
    add_dependencies(all_executables sage_executable)

    # Create a custom target to install sage executable
    add_custom_target(install_sage
        COMMAND ${CMAKE_COMMAND} -E copy_if_different 
                ${CMAKE_SOURCE_DIR}/sage-model/sage 
                ${CMAKE_BINARY_DIR}/sage
        DEPENDS sage_executable
        COMMENT "Installing SAGE executable to build directory"
        VERBATIM
    )

    # Clean target for sage
    add_custom_target(clean_sage
        COMMAND make clean
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/sage-model
        COMMENT "Cleaning SAGE build files"
        VERBATIM
    )

    message(STATUS "SAGE executable target added successfully.")
    endif()  # End of Makefile check
else()
    message(WARNING "sage-model directory not found. SAGE executable will not be built.")
    message(WARNING "To enable SAGE build, run: git submodule update --init --recursive")
endif()

# Enable testing
enable_testing()
#include(tests)

# Install the lightcone.sh script to the runtime directory
install(PROGRAMS scripts/lightcone.sh DESTINATION . 
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

# Install the h5_to_csv executable to runtime/bin
#install(TARGETS h5_to_csv DESTINATION bin)

# Print status message
message(STATUS "Finished initializing custom targets.")
message(STATUS "")
