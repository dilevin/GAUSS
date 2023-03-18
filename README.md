<img src="GaussLogo.png" height="150"></img> <br>
OSX (el capitan)/Ubuntu (14.04): [![Build Status](https://travis-ci.org/dilevin/GAUSS.svg?branch=master)](https://travis-ci.org/dilevin/GAUSS) Windows (Visual Studio 2017): [![Build status](https://ci.appveyor.com/api/projects/status/oyvda3s704ibkfer?svg=true)](https://ci.appveyor.com/project/dilevin/gauss)

### External Dependencies  ###
1. Qt >= 5.8 (https://www.qt.io, version 5.9.0 only for Windows build)
2. OPTIONAL: Pardiso Solver (http://www.pardiso-project.org)
3. OPTIONAL: Gurobi Solver (http://www.gurobi.com)

### Included Submodules (Installed Automatically) ###
1. Libigl https://github.com/libigl/libigl
2. Eigen >= 3.2 (By Default GAUSS uses the libigl Eigen install)
3. Spectra Eigenproblem Solver (https://spectralib.org)
4. The Flexible Collision Library (https://github.com/flexible-collision-library/fcl)
5. Eigen-Gurobi Interface (https://github.com/jrl-umi3218/eigen-gurobi/)

### New Sample Application ###
Use this [sample application](https://github.com/dilevin/GaussApplication/) to get a fast start on developping applications using GAUSS. It provides a working (I think ?) CMake file and includes GAUSS as a submodule.

### "Easy" Install Scripts ###
Clone this repository with 
```bash
git clone --recursive https://github.com/dilevin/GAUSS.git
```

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

#### Build (Linux and OSX) ####

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
	./bin/Example4 # If OpenMP was enabled (OSX only)
	
#### Build (Windows) ####

In the root of this repository do the following

    mkdir build
    cd build
    cmake .. -DUSE_MATLAB=True
    cmake --build . --config Release

To see if it worked, run the tests and examples
	
	./bin/Release/Tests.exe
	./bin/Release/Example1.exe
	./bin/Release/Example2.exe
	./bin/Release/Example3.exe

#### A Note About the Examples ####
To play examples in Gauss press 'p' once the Qt window appears.

#### MATLAB Interface ####
GAUSS Includes a rudimentary MATLAB interface, tested using MATLAB 2015b and 2017a on OSX. To enable the MATLAB interface, build the project Gauss_MATLAB, then open MATLAB and issue the following commands:

	addpath('{Gauss_Root_Dir}/src/MATLAB/')
	addpath('{Gauss_Root_Dir}/build/lib/{Build_Mode_of_Gauss_MATLAB}/')
	savepath

On Windows, an additional path needs to be added:

	addpath('{Gauss_Root_Dir}/build/bin/{Build_Mode_of_Gauss_MATLAB}/')

An example of using the MATLAB interface is given in `{Gauss_Root_Dir}/src/Examples/example8.m`
