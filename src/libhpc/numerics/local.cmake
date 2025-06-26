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
list(APPEND INCFILES "clip.hh" )
list(APPEND INCFILES "constants.hh" )
list(APPEND INCFILES "coords.hh" )
list(APPEND INCFILES "dft.hh" )
list(APPEND INCFILES "eigen_matrix.hh" )
list(APPEND INCFILES "eigen_vector.hh" )
list(APPEND INCFILES "generators.hh" )
list(APPEND INCFILES "histogram.hh" )
list(APPEND INCFILES "integrate.hh" )
list(APPEND INCFILES "interp.hh" )
list(APPEND INCFILES "interp_iterator.hh" )
list(APPEND INCFILES "iqr.hh" )
list(APPEND INCFILES "kde.hh" )
list(APPEND INCFILES "least_squares.hh" )
list(APPEND INCFILES "matrix.hh" )
list(APPEND INCFILES "polynomial.hh" )
list(APPEND INCFILES "quadrature.hh" )
list(APPEND INCFILES "simpson.hh" )
list(APPEND INCFILES "spline.hh" )
list(APPEND INCFILES "spline_integrator.hh" )
list(APPEND INCFILES "tridiag.hh" )
list(APPEND INCFILES "vector.hh" )

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
