# NOTE: This file is marked as assume-unchanged so that users may set custom paths without commiting to 
# the repo. If you would like to make a permanent change to this file, please use the commands
# $ git update-index --no-assume-unchanged config.cmake 
# <Make and commit your changes>
# $ git update-index --assume-unchanged config.cmake


### libigl
set(LIBIGL_INCLUDE_PATH "../libigl/include" CACHE FILEPATH "Root include directory for libigl")


### OpenMP config
set(USE_OPENMP OFF CACHE BOOL "OpenMP requires use of llvm via homebrew on OSX")

# These need only be set/updated if you are trying to use OpenMP on OSX and have installed llvm via homebrew
set(LLVM_BIN "/usr/local/opt/llvm/bin" CACHE STRING "CLANG Binary Directory such as /usr/local/opt/llvm/bin")
set(LLVM_LIB "/usr/local/opt/llvm/lib" CACHE STRING "CLANG Lib Directory such as /usr/local/opt/llvm/lib")
set(LLVM_INCLUDE "/usr/local/opt/llvm/lib/clang/4.0.1/include" CACHE STRING "CLANG Header Directory such as /usr/local/opt/llvm/include")


### Paradiso
set(USE_PARDISO OFF CACHE BOOL "Use Pardiso if it is available to you")

# If USE_PARDISO is ON then we need a path to find it
set(PARDISO_LIB "/usr/local/lib/libpardiso500-MACOS-X86-64.dylib" CACHE STRING "Pardiso Library to use")
