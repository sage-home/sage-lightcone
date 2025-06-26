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
list(APPEND INCFILES "anon.hh" )
list(APPEND INCFILES "application.hh" )
list(APPEND INCFILES "array.hh" )
list(APPEND INCFILES "assign.hh" )
list(APPEND INCFILES "cc_version.hh" )
list(APPEND INCFILES "cuda.hh" )
list(APPEND INCFILES "deallocate.hh" )
list(APPEND INCFILES "file_descriptor.hh" )
list(APPEND INCFILES "filesystem.hh" )
list(APPEND INCFILES "functional.hh" )
list(APPEND INCFILES "has.hh" )
list(APPEND INCFILES "id.hh" )
list(APPEND INCFILES "main.hh" )
list(APPEND INCFILES "mask.hh" )
list(APPEND INCFILES "math.hh" )
list(APPEND INCFILES "matrix.hh" )
list(APPEND INCFILES "narg.hh" )
list(APPEND INCFILES "path_finder.hh" )
list(APPEND INCFILES "random.hh" )
list(APPEND INCFILES "reallocate.hh" )
list(APPEND INCFILES "stream.hh" )
list(APPEND INCFILES "string.hh" )
list(APPEND INCFILES "timer.hh" )
list(APPEND INCFILES "timer_handle.hh" )
list(APPEND INCFILES "tmpfile.hh" )
list(APPEND INCFILES "type_traits.hh" )
list(APPEND INCFILES "varray.hh" )
list(APPEND INCFILES "view.hh" )

# Add source files
# eg. list(APPEND SRCFILES "file" )

# **********************************************************
# ** The template places all source files here, but some/all
# ** may need to be manually placed in the EXEFILES list
list(APPEND SRCFILES "application.cc" )
list(APPEND SRCFILES "file_descriptor.cc" )
list(APPEND SRCFILES "filesystem.cc" )
list(APPEND SRCFILES "id.cc" )
list(APPEND SRCFILES "path_finder.cc" )
list(APPEND SRCFILES "random.cc" )
list(APPEND SRCFILES "stream.cc" )
list(APPEND SRCFILES "tmpfile.cc" )
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
