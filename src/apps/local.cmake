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
#list(APPEND SRCDIRS "analytic" )
#list(APPEND SRCDIRS "bolshoi" )
#list(APPEND SRCDIRS "cli" )
#list(APPEND SRCDIRS "dbcheck" )
#list(APPEND SRCDIRS "dust_plots" )
#list(APPEND SRCDIRS "infoh5" )
#list(APPEND SRCDIRS "magnitudes" )
#list(APPEND SRCDIRS "rebin" )
#list(APPEND SRCDIRS "ssp2h5" )
#list(APPEND SRCDIRS "ssp_restrict" )
#list(APPEND SRCDIRS "subcones" )
#list(APPEND SRCDIRS "validate" )
#list(APPEND SRCDIRS "zen" )
# **********************************************************

# Add header files
# eg. list(APPEND INCFILES "file" )
list(APPEND INCFILES "kdapplication.hh" )

# Add source files
# eg. list(APPEND SRCFILES "file" )

# **********************************************************
# ** The template places all source files here, but some/all
# ** may need to be manually placed in the EXEFILES list
#list(APPEND SRCFILES "add_tree_info.cc" )
list(APPEND SRCFILES "kdapplication.cc" )
# **********************************************************

# Add executable source files (those with a main())
# eg. list(APPEND EXEFILES "file" )
list(APPEND EXEFILES "cli_lightcone.cc" )
list(APPEND EXEFILES "sage2h5.cc" )
list(APPEND EXEFILES "sage2kdtree.cc" )
list(APPEND EXEFILES "dstreeinit.cc" )
list(APPEND EXEFILES "sageh5tokdtree.cc" )
list(APPEND EXEFILES "sageh5toxml.cc" )
list(APPEND EXEFILES "sageh5toh5.cc" )
list(APPEND EXEFILES "sageimport.cc" )

# Add data files
# eg. list(APPEND DATAFILES "file" )

# Set data subdirectory
# eg. set(DATASUBDIR "dir" )

# Add subdirectories that we want to decend into
#   but which we won't scan for sources, etc
# eg. list(APPEND PASSDIRS "dir" )
