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
# **********************************************************

# Add header files
# eg. list(APPEND INCFILES "file" )
list(APPEND INCFILES "application.hh" )
list(APPEND INCFILES "assert.hh" )
list(APPEND INCFILES "async.hh" )
list(APPEND INCFILES "comm.hh" )
list(APPEND INCFILES "comm_proxies.hh" )
list(APPEND INCFILES "cuda.hh" )
list(APPEND INCFILES "datatype.hh" )
list(APPEND INCFILES "helpers.hh" )
list(APPEND INCFILES "host.hh" )
list(APPEND INCFILES "indexer.hh" )
list(APPEND INCFILES "init.hh" )
list(APPEND INCFILES "insist.hh" )
list(APPEND INCFILES "logger.hh" )
list(APPEND INCFILES "main.hh" )
list(APPEND INCFILES "managed.hh" )
list(APPEND INCFILES "partition.hh" )
list(APPEND INCFILES "request.hh" )
list(APPEND INCFILES "requests.hh" )
list(APPEND INCFILES "scatter.hh" )
list(APPEND INCFILES "surjection.hh" )
list(APPEND INCFILES "type_map.hh" )
list(APPEND INCFILES "vct.hh" )

# Add source files
# eg. list(APPEND SRCFILES "file" )

# **********************************************************
# ** The template places all source files here, but some/all
# ** may need to be manually placed in the EXEFILES list
list(APPEND SRCFILES "application.cc" )
list(APPEND SRCFILES "async.cc" )
list(APPEND SRCFILES "comm.cc" )
list(APPEND SRCFILES "cuda.cc" )
list(APPEND SRCFILES "datatype.cc" )
list(APPEND SRCFILES "helpers.cc" )
list(APPEND SRCFILES "host.cc" )
list(APPEND SRCFILES "init.cc" )
list(APPEND SRCFILES "logger.cc" )
list(APPEND SRCFILES "partition.cc" )
list(APPEND SRCFILES "request.cc" )
list(APPEND SRCFILES "requests.cc" )
list(APPEND SRCFILES "scatter.cc" )
list(APPEND SRCFILES "surjection.cc" )
list(APPEND SRCFILES "type_map.cc" )
list(APPEND SRCFILES "vct.cc" )
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
