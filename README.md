# ASVO-TAO Science Modules #

## Description

These are the science modules, designed to carry out the HPC operations, such
as building a lightcone and calculating the SEDs.

## Dependencies

 * Python - An object-oriented interpreted programming language.
 * Boost - A set of efficient C++ extensions.
 * MPI - Message Passing Interface implementation. Tested with OpenMPI.
 * HDF5 - Hierarchical Data Format libraries. Currently tested with version 1.8.7.
    (Must be compiled with OpenMPI or compilation errors will be encountered)
 * GSL - The GNU Scientific Library.
 * PugiXML - A C++ XML parsing library.
 * cfitsio - A C implementation of a FITS file reader/writer.
 * SOCI - A C++ database access wrapper, capable of using a wide variety of drivers.
 * PostgresQL (optional) - A production quality database management system.
 * sqlite3 (optional) - A simple SQL database implementation, used for testing only.
 
 * git submodules - You may need to run `git submodule update --init` after a git clone to get these.

## Building

The science modules use autotools as the build system. As with any package
that has numerous dependencies, we recommend installing as many of them
as possible using your system's package management software.
Note that some of the dependencies are included as git submodules.
Ensure you run git submodule update to check them out.

Once you have your dependencies in place, run

```bash
cmake -DCMAKE_BUILD_TYPE=Debug . && make
```

to build the science modules. If you find that there are some dependencies
that could not be found you can specify them on the command line to `./configure` using
various options. For example, to specify the location of Boost you may
give the base location:

```bash
cmake -DHDF5_ROOT=/usr/lib/x86_64-linux-gnu/hdf5/openmpi
```

Note that you can also supply multiple directories by separating
them with a colon:

```bash
./configure --with-boost_inc=/path/to/boost=/path1:/path2 
```

To see more options, please run:

```bash
./configure --help
```

## Testing

To run the full suite of unit-tests, please run:
... TODO (unsure how to run these post autotools, previous it was `scons check`)

## Installing

To install, please use:

```bash
make install
```

You may specify the installation path using `PREFIX` when configuring:

```bash
./configure PREFIX=/path/to/install 
make
make install
```

## G2 installation
A compilation script is provided for the g2 environment which loads the required modules, runs `bootstrap` `configure` and `make`.

This script may be useful as a reference for how to use some of the `./configure` options as well.

```bash
./sstar-compile.sh`
```

## Development Environment (docker)
There are a few docker images that exist which can make local development much easier.
As a vague reference you can look at the CircleCI configuration (which is not passing due to some failing unit tests, but the docker instance does work) in `.circle/config.yml`

### Set up:
* `git clone https://github.com/ASVO-TAO/TAO-ScienceModules`
* `git submodule update --init`
* `pip install -qr requirements.txt`
* `fab build`

### To run unit tests
* `ctest`

### To run a specific job
These will require preparing an XML file for the job (The `params.xml` file that exists in the job directory on TAO, note this is *NOT* the same as the `params.xml` file that the UI lets you download.) Ensure that the `<database>` key in your params matches the name of the database on your docker db. By default the docker db has minimillennium is installed in the `postgres` database.

You may also need to set up a database for this job if you're not using mini-millennium (which docker has by default). You can connect to the database with `fab psql` if you need to view it. (Or `fab cli` to get a bash prompt then `psql -h db -U postgres postgres`)

Note that you cannot be running a postgresql server on the docker host, or the ports will conflict.
* `fab run:<params.xml>` (eg. `fab run:params.xml`)
* To use the regression testing xml file as an example `fab run:./regression/unique_cones/params.xml`

To debug it
* `fab gdb:<params.xml>` or `fab gdb_remote:<params.xml>`

The latter can be useful when using an IDE such as cLion (you can run a remote GDB session in cLion which docker will then connect to)

`fabfile.py` has various other useful functions for other things you may want to do.


### Technical note on format of skymaker ".list" file
Bulge:
200, px, py, mag, bulge_to_total,                                         scale_radius, 1.0, pos_angle,          0.0, 1.0, pos_angle, redshift
Disk:
200, px, py, mag, bulge_to_total,                                                  0.0, 1.0, pos_angle, scale_radius,  ar, pos_angle, redshift
Disk_Bulge:
200, px, py, mag, bulge_to_total, get_bulge_to_disk_scale(bulge_to_total)*scale_radius, 1.0, pos_angle, scale_radius,  ar, pos_angle, redshift

