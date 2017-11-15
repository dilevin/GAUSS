[![Build Status](https://travis-ci.com/dilevin/GAUSS.svg?token=poS417w5DgfGAsnmYggm&branch=master)](https://travis-ci.com/dilevin/GAUSS)

### Generic Build Instructions ###

#### Setup ####

Ensure that Qt version >= 5.8 is installed.

Download (and compile?) the igl library into the directory above `${GAUSS_DIR}` or choose a custom path in `config.cmake`.

If you are on OSX and you wish to enable OpenMP, install llvm via homebrew and follow the instructions in `config.cmake`

Likewise, if you wish to use the pardiso solver, download the pardiso library from the website and set the path in `config.cmake`. If you get errors referring to libgfortran, you may need to follow [this](http://www.alecjacobson.com/weblog/?p=3946) guide to fix the library for your paths.

#### Build ####

In the root of this repository do the following

    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release -C ../config.cmake ../
    make

To see if it worked, run the tests and examples
	
	./bin/Tests
	./bin/Example1
	./bin/Example2
	./bin/Example3
	./bin/Example4 # If OpenMP was enabled

#### MATLAB Interface ####
GAUSS Includes a rudimentary MATLAB interface, tested using MATLAB 2015b on OSX. To enable the MATLAB interface, build the project Gauss_MATLAB, then open MATLAB and issue the following commands:
	addpath('{Gauss_Root_Dir}/src/MATLAB/')
	addpath('{Gauss_Root_Dir}/build/lib/{Build_Mode_of_Gauss_MATLAB}/')
	savepath

An example of using the MATLAB interface is given in {Gauss_Root_Dir}/src/Examples/example8.m
