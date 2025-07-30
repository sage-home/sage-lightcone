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
list(APPEND SRCDIRS "kdtree_backend" )
list(APPEND SRCDIRS "utils" )
# **********************************************************

# Add header files
# eg. list(APPEND INCFILES "file" )
list(APPEND INCFILES "backend.hh" )
list(APPEND INCFILES "base.hh" )
list(APPEND INCFILES "batch.hh" )
list(APPEND INCFILES "box.hh" )
list(APPEND INCFILES "factory.hh" )
list(APPEND INCFILES "filter.hh" )
list(APPEND INCFILES "galaxy_iterator.hh" )
list(APPEND INCFILES "globals.hh" )
list(APPEND INCFILES "integration.hh" )
list(APPEND INCFILES "lightcone.hh" )
list(APPEND INCFILES "lightcone_galaxy_iterator.hh" )
list(APPEND INCFILES "lightcone_tile_iterator.hh" )
list(APPEND INCFILES "module.hh" )
list(APPEND INCFILES "query.hh" )
list(APPEND INCFILES "simulation.hh" )
list(APPEND INCFILES "subcones.hh" )
list(APPEND INCFILES "tile.hh" )
list(APPEND INCFILES "types.hh" )
list(APPEND INCFILES "utils.hh" )
list(APPEND INCFILES "xml_dict.hh" )
list(APPEND INCFILES "data_dict.hh" )

# Add source files
# eg. list(APPEND SRCFILES "file" )

# **********************************************************
# ** The template places all source files here, but some/all
# ** may need to be manually placed in the EXEFILES list
list(APPEND SRCFILES "backend.cc" )
list(APPEND SRCFILES "filter.cc" )
list(APPEND SRCFILES "globals.cc" )
list(APPEND SRCFILES "lightcone.cc" )
list(APPEND SRCFILES "lightcone_tile_iterator.cc" )
list(APPEND SRCFILES "simulation.cc" )
list(APPEND SRCFILES "subcones.cc" )
list(APPEND SRCFILES "utils.cc" )
list(APPEND SRCFILES "xml_dict.cc" )
list(APPEND SRCFILES "data_dict.cc" )
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
