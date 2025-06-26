# Configure any external projects, submodules, etc.
macro(local_extern cur_dir )
    # Include external libraries added as submodules to the build here, by
    # using 'add_external_submodule' like follows:
    #
    #    add_external_submodule( ${cur_dir} "library_name" ${cur_dir}"/library_name/file_to_check" )
    #
    # where file to check is a file that will always be in the project (eg. README.md)
endmacro()
