# *****************************************************
# *** A 'local.cmake' file following this format    ***
# *** should reside in every directory contributing ***
# *** to the C/C++ code of this project.            ***
# *****************************************************

# The following directive is used to organise
# how the Doxygen groups defined by this project
# will be organised in the documentation:
#
# set_active_API_module root Project root

# Set empty defaults
set(LIBDIRS    "" )
set(SRCDIRS    "" )
set(INCFILES   "" )
set(SRCFILES   "" )
set(EXEFILES   "" )
set(DATAFILES  "" )
set(PASSDIRS   "" )
set(DATASUBDIR "" )

# Add subdirectories that are roots to libraries
# eg. list(APPEND LIBDIRS "dir" )

# Add directories that contribute source files 
# eg. list(APPEND SRCDIRS "dir" )

# **********************************************************
# ** The template places all directories here, but some/all
# ** may need to be manually placed in the LIBDIRS list
list(APPEND SRCDIRS "algorithm" )
list(APPEND SRCDIRS "containers" )
list(APPEND SRCDIRS "debug" ) # Files are missing from this path including 'omp_help.hh'
list(APPEND SRCDIRS "h5" )
list(APPEND SRCDIRS "logging" )
list(APPEND SRCDIRS "mpi" )
list(APPEND SRCDIRS "numerics" )
list(APPEND SRCDIRS "system" ) # Files are missing from this path including 'epoll.h'
# **********************************************************

# Add header files
# eg. list(APPEND INCFILES "file" )
list(APPEND INCFILES "algorithm.hh" )
list(APPEND INCFILES "debug.hh" )
list(APPEND INCFILES "h5.hh" )
list(APPEND INCFILES "libhpc.hh" )
list(APPEND INCFILES "logging.hh" )
list(APPEND INCFILES "mpi.hh" )
list(APPEND INCFILES "numerics.hh" )
list(APPEND INCFILES "system.hh" )

# Add source files
# eg. list(APPEND SRCFILES "file" )

# **********************************************************
# ** The template places all source files here, but some/all
# ** may need to be manually placed in the EXEFILES list
# **********************************************************

# Add executable source files (those with a main())
# eg. list(APPEND EXEFILES "file" )

# Add data files
# eg. list(APPEND DATAFILES "file" )

# Set data subdirectory
# eg. set(DATASUBDIR "dir" )

# Add subdirectories that we want to decend into
#   but which we won't scan for sources, etc
# eg. list(APPEND PASSDIRS "dir" )
