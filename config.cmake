
### libigl
set(LIBIGL_INCLUDE_PATH "/Users/lawson/Workspace/dgp/libigl/include" CACHE FILEPATH "Root include directory for libigl")


### OpenMP config
set(USE_OPENMP OFF CACHE BOOL "OpenMP requires use of llvm via homebrew on OSX")

# These need only be set if you are trying to use OpenMP on OSX 
set(LLVM_BIN "/usr/local/opt/llvm/bin" CACHE STRING "CLANG Binary Directory")
set(LLVM_LIB "/usr/local/opt/llvm/lib" CACHE STRING "CLANG Lib Directory")
set(LLVM_INCLUDE "/usr/local/opt/llvm/include" CACHE STRING "CLANG Header Directory")