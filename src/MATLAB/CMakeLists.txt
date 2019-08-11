cmake_minimum_required(VERSION 3.2)

project(Gauss_MATLAB)

find_package(MATLAB REQUIRED)

include_directories(${Matlab_INCLUDE_DIRS})

get_filename_component(MATLAB_PATH_DIR "${Matlab_MEX_LIBRARY}" PATH)

#have to manually ad libmx on OSX (on my machine at least) for some reason
find_library(LUA_LIB lua)
if(NOT Matlab_MX_LIBRARY)
    message(WARNING "MX Library not found, setting by hand")
    if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
        set(Matlab_MX_LIBRARY ${MATLAB_PATH_DIR}/libmx.lib)
    elseif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
        set(Matlab_MX_LIBRARY ${MATLAB_PATH_DIR}/libmx.so)
    elseif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        if(APPLE)
            set(Matlab_MX_LIBRARY ${MATLAB_PATH_DIR}/libmx.dylib)
        endif(APPLE)
    endif()

endif()

if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
     matlab_add_mex(
        NAME Gauss_MATLAB
        SRC class_handle.hpp WorldFEM_Interface.cpp UtilitiesEigenMex.h
        LINK_TO ${Matlab_MX_LIBRARY} ${Matlab_MEX_LIBRARY} Core Base SolversLinear ParticleSystem FEM ${Matlab_ROOT_DIR}/sys/os/glnxa64/libiomp5.so
        OUTPUT_NAME WorldFEM_Interface)
elseif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	if(APPLE)
		matlab_add_mex(
			NAME Gauss_MATLAB 
			SRC class_handle.hpp WorldFEM_Interface.cpp UtilitiesEigenMex.h
			LINK_TO ${Matlab_MX_LIBRARY} ${Matlab_MEX_LIBRARY} Core Base SolversLinear ParticleSystem FEM  ${Matlab_ROOT_DIR}/sys/os/maci64/libiomp5.dylib
			OUTPUT_NAME WorldFEM_Interface)
	endif(APPLE)
elseif (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    matlab_add_mex(
        NAME Gauss_MATLAB 
        SRC class_handle.hpp WorldFEM_Interface.cpp UtilitiesEigenMex.h
    LINK_TO ${Matlab_MX_LIBRARY} ${Matlab_MEX_LIBRARY} Core Base SolversLinear ParticleSystem FEM  ${Matlab_ROOT_DIR}/bin/win64/libiomp5md.lib
        OUTPUT_NAME WorldFEM_Interface)
endif()

