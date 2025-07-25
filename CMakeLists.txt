# This is the CMakeList.txt file used to build the project
#    - We use this old version because it is the  -
#    -   latest version availale on Travis CI     -
cmake_minimum_required(VERSION 3.2)

# Read the file with the macro definitions which specify
# the behaviour specific to this project and then run
# the initialization macro
include(${CMAKE_SOURCE_DIR}/project.cmake)

# Add the 'cmake' directory in the project directory to the module path
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_SOURCE_DIR};${CMAKE_SOURCE_DIR}/.cmake")

# Load and run a macro to validate the initialization
include ("init_validate")
init_validate()

# Set install directory to ./runtime
# (user can still override with DESTDIR with make to set installation directory)
set( CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/runtime" )

# Set default compilers (must come before project line)
# For GCC build
set(CMAKE_C_COMPILER   "gcc")
set(CMAKE_CXX_COMPILER "g++")

# For Clang build
#set(CMAKE_C_COMPILER   "clang")
#set(CMAKE_CXX_COMPILER "clang++")

# Compiler flags - separate C and C++ flags
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu99 -Wextra -mtune=native -march=native")
#set(CMAKE_CXX_FLAGS "-std=c++14 -Wextra -mtune=native -march=native")
#set(CMAKE_CXX_FLAGS "-std=c++17 -Wextra -mtune=native -march=native")
#set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_EXTENSIONS OFF)
# Configure the cmake C++ compile features
# To change the C++ standard at build time (to c++23, for example), run cmake like so:
#    `cmake -DCMAKE_CXX_STANDARD=23 .."
set(CXX_STANDARD_DEFAULT 20)
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
add_definitions("-DPARALLELHDF5")

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
else()
    message(STATUS "Compiling in RELEASE mode...")
    add_definitions(-DNDEBUG -DNLOGDEBUG -DNLOGTRIVIAL -DNINSTRUMENT -DNSTACKTRACE -DNMEMDEBUG -DNMEMOPS -DNMEMSTATS)

    # CURRENTLY NDEBUG BREAKS THE PROJECT AT RUNTIME
    # TODO: Remove this once we find out why NDEBUG breaks at runtime
    #message(FATAL_ERROR "COMPILING IN RELEASE MODE IS KNOWN TO BREAK AT RUNTIME. -DNDEBUG IS BROKEN.")

    # Release compiler flags
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O1")
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
