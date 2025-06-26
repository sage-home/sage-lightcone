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
#list(APPEND SRCDIRS "thrust" )
# **********************************************************

# Add header files
# eg. list(APPEND INCFILES "file" )
list(APPEND INCFILES "bin.hh" )
list(APPEND INCFILES "binary_partitioner.hh" )
list(APPEND INCFILES "brent.hh" )
list(APPEND INCFILES "counts.hh" )
list(APPEND INCFILES "cuda_helpers.hh" )
list(APPEND INCFILES "dual.hh" )
list(APPEND INCFILES "inner_product.hh" )
list(APPEND INCFILES "kdtree.hh" )
list(APPEND INCFILES "median.hh" )
list(APPEND INCFILES "morton.hh" )
list(APPEND INCFILES "multimatch.hh" )
list(APPEND INCFILES "newton.hh" )
list(APPEND INCFILES "permute.hh" )
list(APPEND INCFILES "reorder.hh" )
list(APPEND INCFILES "richardson.hh" )
list(APPEND INCFILES "ridders.hh" )
list(APPEND INCFILES "select.hh" )
list(APPEND INCFILES "sequence.hh" )
list(APPEND INCFILES "sequence_checks.hh" )
list(APPEND INCFILES "sort_by_key.hh" )
list(APPEND INCFILES "sort_permute2_iter.hh" )
list(APPEND INCFILES "sort_permute_iter.hh" )
list(APPEND INCFILES "split.hh" )
list(APPEND INCFILES "standard.hh" )
list(APPEND INCFILES "uniform.hh" )
list(APPEND INCFILES "xdmf_writer.hh" )

# Add source files
# eg. list(APPEND SRCFILES "file" )

# **********************************************************
# ** The template places all source files here, but some/all
# ** may need to be manually placed in the EXEFILES list
list(APPEND SRCFILES "morton.cc" )
list(APPEND SRCFILES "multimatch.cc" )
list(APPEND SRCFILES "xdmf_writer.cc" )
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
