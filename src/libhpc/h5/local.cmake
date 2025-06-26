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
list(APPEND INCFILES "buffer.hh" )
list(APPEND INCFILES "dataset.hh" )
list(APPEND INCFILES "dataspace.hh" )
list(APPEND INCFILES "datatype.hh" )
list(APPEND INCFILES "derive.hh" )
list(APPEND INCFILES "file.hh" )
list(APPEND INCFILES "group.hh" )
list(APPEND INCFILES "location.hh" )
list(APPEND INCFILES "property_list.hh" )
list(APPEND INCFILES "types.hh" )

# Add source files
# eg. list(APPEND SRCFILES "file" )

# **********************************************************
# ** The template places all source files here, but some/all
# ** may need to be manually placed in the EXEFILES list
list(APPEND SRCFILES "buffer.cc" )
list(APPEND SRCFILES "dataset.cc" )
list(APPEND SRCFILES "dataspace.cc" )
list(APPEND SRCFILES "datatype.cc" )
list(APPEND SRCFILES "derive.cc" )
list(APPEND SRCFILES "file.cc" )
list(APPEND SRCFILES "group.cc" )
list(APPEND SRCFILES "location.cc" )
list(APPEND SRCFILES "property_list.cc" )
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
