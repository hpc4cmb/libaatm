
# libaatm

[![Build Status](https://travis-ci.org/hpc4cmb/libaatm.svg?branch=master)](https://travis-ci.org/hpc4cmb/libaatm)

This software is a repackaged version of the Alma Atmospheric Transmission at
Microwaves Tools.  The original upstream code came from this git repository:

https://bitbucket.sco.alma.cl/projects/ASW/repos/telcal/browse/TelCalResults/Libraries/ATM

The git version of that source is recorded in the VERSION text file in this
current directory.

Only minor changes were made to source files and tests to resolve compiler
warnings, etc.

## Installation

This library uses a cmake build system.  You should make a build directory:

    %>  mkdir build
    %>  cd build

And then run cmake with the options you wish (like the install location):

    %>  cmake -DCMAKE_INSTALL_PREFIX=/path/to/install ..

Next build with:

    %>  make -j 4

Run tests with:

    %>  make test

And install the package with:

    %>  make install
