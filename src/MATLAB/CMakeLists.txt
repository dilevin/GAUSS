cmake_minimum_required(VERSION 3.5)

project(Gauss_MATLAB)

find_package(MATLAB REQUIRED)

include_directories(${Matlab_INCLUDE_DIRS})

message(WARNING ${Matlab_LIBRARIES})
message(WARNING ${Matlab_MX_LIBRARY})

get_filename_component(MATLAB_PATH_DIR "${Matlab_MEX_LIBRARY}" PATH)

#have to manually ad libmx on OSX (on my machine at least) for some reason
matlab_add_mex(
    NAME Gauss_MATLAB 
    SRC class_handle.hpp WorldFEM_Interface.cpp UtilitiesEigenMex.h
    LINK_TO ${MATLAB_PATH_DIR}/libmx.dylib ${Matlab_MEX_LIBRARY} Core Base SolversLinear ParticleSystem FEM  ${LLVM_LIB}/libiomp5.dylib
    OUTPUT_NAME WorldFEM_Interface
)

