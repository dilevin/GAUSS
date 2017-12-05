<img src="GaussLogo.png" height="150"></img> <br>
OSX (el capitan)/Ubuntu (14.04): [![Build Status](https://travis-ci.org/dilevin/GAUSS.svg?branch=master)](https://travis-ci.org/dilevin/GAUSS) Windows (Visual Studio 2017): [![Build status](https://ci.appveyor.com/api/projects/status/oyvda3s704ibkfer?svg=true)](https://ci.appveyor.com/project/dilevin/gauss)

### Dependencies ###
1. Eigen >= 3.2
2. Libigl https://github.com/libigl/libigl
3. Qt >= 5.8 (https://www.qt.io, version 5.9.0 only for Windows build)
4. OPTIONAL: Pardiso Solver (http://www.pardiso-project.org)
5. OPTIONAL: Spectra Eigenproblem Solver (https://spectralib.org)

### "Easy" Install Scripts ###
#### Install Instructions OS X (Basic install, no bells or whistles) ####
1. Install Homebrew (https://brew.sh)
2. At command line: chmod a+x ./InstallGAUSS_OSX.sh
3. Run InstallGAUSS_OSX.sh which does the following:
	- downloads and installs Qt 5.9 in ~/Qt
	- installs CMake using homebrew (upgrades currently installed version)
	- installs Eigen 3 using home brew (upgrades currently installed version)
	- installs libigl in ./build/libigl
	- installs Gauss in ./build

#### Install Instructions Ubuntu (Basic install, no bells or whistles) ####
1. This install procedure requires gcc and g++ version 5 or greater setup as the default c/c++ compilers
2. At command prompt: chmod a+x ./InstallGAUSS_Ubuntu.sh
3. Run InstallGAUSS_Ubuntu.sh which does the following
	- downloads and installs Qt 5.9 in ~/Qt
	- install libigl into ./build/libigl
	- installs CMake using apt-get
	- installs Eigen 3 using apt-get
	- installs Gauss in ./build

#### Warning ####
These scripts are hand maintained to make GAUSS Easy to build on certain platforms. If they fail please resort to custom install using cmake below. 

### Generic Build Instructions ###

#### Setup ####

Ensure that Qt version >= 5.8 is installed.
On Ubuntu 16.04 you may have to change the paths to be consistent in ccmake. Something like the following.
    /home/<user>/Qt/5.9.2/gcc_64/lib/cmake/

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

#### Notes for Visual Studio ####
To build the MATLAB interface using Visual Studio, currnetly you have to manually modify a couple of linker inputs. After building the VS solution using CMake, open it in VS, and go to the linker input for the Gauss_MATLAB project. That is, Gauss_MATLAB (Solution Explorer) > Properties > Linker > Input > Additional Dependencies.
Change `libmx.dylib` to `libmx.lib`, and remove the input `libiomp5.lib`. Now you should be able to build hte MATLAB interface on Windows.
