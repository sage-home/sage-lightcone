# TAO LIGHTCONE CLI #

## Description

These are the science modules, designed to carry out the HPC operation of extracting lightcones from "KD TREE" indexed simulations.

## Dependencies

 * Boost - A set of efficient C++ extensions.
 * MPI - Message Passing Interface implementation. Tested with OpenMPI.
 * HDF5 - Hierarchical Data Format libraries. Currently tested with version 1.8.7.
    (Must be compiled with OpenMPI or compilation errors will be encountered)
 * GSL - The GNU Scientific Library.
 * PugiXML - A C++ XML parsing library.
 

## Building

To build on Swinburne's nt.swin.edu.au after cloning from the repository just run
```bash
source ./setup.sh
./run_build.sh
```

To build in other environments refer to the .gitlab-ci.yml on how this is done within a CI/CDD framework. The science modules use cmake as a build system. As with any package
that has numerous dependencies, we recommend installing as many of them
as possible using your system's package management software. Note that the source to a compatible puxixml library is already included under the dep directory.

## Testing

The script test_cli_validation.sh will test that the build has succeeded and that a suite of validation tests correctly inform the user of errors in the command line arguments.

## Installing

```bash
cd bin
make install
```

The script lightcone.sh is installed under the runtime directory and allows a user to easily run the cli-lightcone tool on nt.swin.edu.au (for members of the oz114 project) ensuring that the right modules are loaded first.


