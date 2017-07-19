# NOTE: This file is marked as assume-unchanged so that users may set custom paths without commiting to 
# the repo. If you would like to make a permanent change to this file, please use the commands
# $ git update-index --no-assume-unchanged config.cmake 
# <Make and commit your changes>
# $ git update-index --assume-unchanged config.cmake


### libigl
set(LIBIGL_INCLUDE_PATH "../libigl/include" CACHE FILEPATH "Root include directory for libigl")

### OpenMP config
set(USE_OPENMP OFF CACHE BOOL "OpenMP requires use of llvm via homebrew on OSX")

# These need only be set if you are trying to use OpenMP on OSX and have installed llvm via homebrew
set(LLVM_BIN "<path>" CACHE STRING "CLANG Binary Directory such as /usr/local/opt/llvm/bin")
set(LLVM_LIB "<path>" CACHE STRING "CLANG Lib Directory such as /usr/local/opt/llvm/lib")
set(LLVM_INCLUDE "<path>" CACHE STRING "CLANG Header Directory such as /usr/local/opt/llvm/include")
