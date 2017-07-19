### Generic Build Instructions ###

#### Setup ####

Ensure that Qt version >= 5.8 is installed.

Download and compile? the igl library into the directory above `${GAUSS_DIR}` or choose a custom path in `config.cmake`.

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
