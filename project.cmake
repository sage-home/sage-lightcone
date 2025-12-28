# Define these outside any macros or functions
set(GBP_PROJECT_NAME  "TAO_ScienceModules")

# Set options and dependencies
macro(project_options_and_dependencies)

    # A bunch of #include paths are expressed relative to 'src', so add it to the list of paths
    list(APPEND INC_DIRS_EXTRA "src" )

    #========= Set project-specific options ========
    # eg. option(SID_DEBUG "Enable SID debugging information" OFF)
    # eg. option(USE_DOUBLE "Use double liberally" OFF)
    #     if(USE_DOUBE)
    #         add_definitions(-DUSE_DOUBLE)
    #     endif()
    
    #=========== Add 3rd-party libraries ===========
    # (look in gbpBuild/cmake/3rd_party.cmake for a list of supported libraries)
    message(STATUS "")
    message(STATUS "Initializing 3rd-party libraries...")
    
    # Three macros have been set-up for automating
    # 3rd-party library configuration.  All accept
    # the described inputs, plus optional arguments
    # for libraries with specific support for it:
    
    # 1) *required* 3rd-Party libraries.
    #
    # First argument is library name.
    #
    # The build will fail if any of these fail to be configured.
    #
    # 2) *requested* 3rd-Party libraries (and their defaults).
    #
    # First argument is library name; second is default state (ON or OFF).
    #
    # The build will proceed if these fail to configure.
    #
    # 3) *optional* 3rd-Party libraries (and their defaults).
    #
    # First argument is library name; second is default state (ON or OFF).
    #
    # The build will succeed if any of these are switched off but
    # will fail if they are turned on and fail to configure.

    # System libraries
    set(THREADS_PREFER_PTHREAD_FLAG ON)
    set_3rd_party_required("Threads")

    # Non-system libraries
    # Note: boost_system is header-only in Boost 1.69+ so not listed as a component
    set(BOOST_COMPONENTS "")
    set(SOCI_COMPONENTS "")
    list(APPEND BOOST_COMPONENTS "program_options")
    list(APPEND BOOST_COMPONENTS "regex")
    list(APPEND BOOST_COMPONENTS "filesystem")
    list(APPEND BOOST_COMPONENTS "chrono")
    list(APPEND BOOST_COMPONENTS "unit_test_framework")

    set_3rd_party_required("Boost" ${BOOST_COMPONENTS})
    set_3rd_party_required("MPI")
    set_3rd_party_required("PugiXML")
    set_3rd_party_required("HDF5")
    set_3rd_party_required("GSL")

    # Set-up documentation support by default
    set_3rd_party_requested("GBP_DOCS_BUILD" ON)

    # Print status message
    message(STATUS "Finished initializing 3rd-party libraries.")
endmacro()
