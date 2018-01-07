cmake_minimum_required (VERSION 3.2)
include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

project (Gauss)

set(CMAKE_CXX_STANDARD 14)
#this code uses OpenMP which is supported almost everywhere except for OS X clang
#to get around this you can choose to install clang via homebrew to get access to OpenMP

#expose USE_OPENMP variable for OS X only
if(APPLE)
    if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        set(USE_OPENMP CACHE BOOL "Use OpenMP via homebrew llvm on OSX")
    endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
endif(APPLE)

#do you want to use Pardiso
set(USE_PARDISO CACHE BOOL "Use Pardiso if available on system")
set(USE_SPECTRA CACHE BOOL "Use Spectra to solve eigenvalue problems (if available)")
set(USE_MATLAB CACHE BOOL "Use OpenMP via homebrew llvm on OSX")

if(APPLE)
    if(USE_OPENMP)
        set(LLVM_BIN CACHE PATH "CLANG Binary Directory")
        set(LLVM_LIB CACHE STRING "CLANG Lib Directory")
        set(LLVM_INCLUDE CACHE STRING "CLANG Header Directory")
    endif(USE_OPENMP)
endif(APPLE)

if(USE_PARDISO)
    set(PARDISO_LIB CACHE STRING "Pardiso Library to use")
endif(USE_PARDISO)

if(USE_SPECTRA)
    set(SPECTRA_INCLUDE_DIR CACHE STRING "Spectra Root Directory")
endif(USE_SPECTRA)

#find_package(OpenMP REQUIRED)
#Global variables to configure other projects
set(GAUSS_INCLUDE_DIR  CACHE INTERNAL "")
set(GAUSS_LINK_DIR  CACHE INTERNAL "")
set(GAUSS_LIBS  CACHE INTERNAL "")
set(GAUSS_BINARY_DIR ${PROJECT_BINARY_DIR} CACHE INTERNAL "")

if(APPLE)
    if(USE_OPENMP)
        set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} -fopenmp)
        set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -fopenmp)
        set(CMAKE_XCODE_ATTRIBUTE_CC /usr/local/opt/llvm/bin/clang)
        set(GAUSS_INCLUDE_DIR ${LLVM_INCLUDE} CACHE INTERNAL "")
        set(GAUSS_LINK_DIR  ${LLVM_LIB} CACHE INTERNAL "")
        set(GAUSS_LIBS ${LLVM_LIB}/libiomp5.dylib CACHE INTERNAL "")
        add_definitions(-DGAUSS_OPENMP)
    endif(USE_OPENMP)
else()
    find_package(OpenMP REQUIRED)
    set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} -fopenmp)
    set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -fopenmp)    
    add_definitions(-DGAUSS_OPENMP)
endif(APPLE)

if(USE_PARDISO)
    set(GAUSS_LIBS ${GAUSS_LIBS} ${PARDISO_LIB} CACHE INTERNAL "")
    get_filename_component(PARDISO_DIR ${PARDISO_LIB} DIRECTORY)
    set(GAUSS_LINK_DIR ${GAUSS_LINK_DIR} ${PARDISO_DIR} CACHE INTERNAL "")
    add_definitions(-DGAUSS_PARDISO)
endif(USE_PARDISO)

if(USE_SPECTRA)
    set(GAUSS_INCLUDE_DIR ${GAUSS_INCLUDE_DIR} ${SPECTRA_INCLUDE_DIR}/include CACHE INTERNAL "")
    add_definitions(-DGAUSS_SPECTRA)
endif(USE_SPECTRA)

#Output directories
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib )

#add define for data directory
add_definitions(-DGAUSS_DATA_DIR=${PROJECT_SOURCE_DIR}/data)

#grab important packages
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

#libigl include directory
set(LIBIGL_INCLUDE_PATH ${PROJECT_SOURCE_DIR}/../libigl/include CACHE FILEPATH "Root include directory for libigl")

#Eigen
find_package(Eigen3 REQUIRED)

#GTest (From Masha Shugrina)
# Google Test (Not using this right now and its causing build problems so comment out)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Download and instalxl
#ExternalProject_Add(
#gtest
#URL https://github.com/google/googletest/archive/release-1.7.0.zip
#PREFIX ${PROJECT_SOURCE_DIR}/ThirdParty/gtest
#INSTALL_DIR ${PROJECT_BINARY_DIR}/lib
#CMAKE_ARGS ${CMAKE_ARGS}
# Disable install step
#INSTALL_COMMAND ""
#                )

# Create a libgtest target to be used as a dependency by test programs
#add_library(libgtest IMPORTED STATIC GLOBAL)
#add_dependencies(libgtest gtest)

# Set gtest properties
#ExternalProject_Get_Property(gtest source_dir binary_dir)
#set_target_properties(libgtest PROPERTIES
#		      "IMPORTED_LOCATION" ${PROJECT_BINARY_DIR}/lib/libgtest.a
#                      "IMPORTED_LINK_INTERFACE_LIBRARIES" "${CMAKE_THREAD_LIBS_INIT}"
#                     )
#Google test ends here

#Individual feature directories
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${LIBIGL_INCLUDE_PATH})
#include_directories(/usr/local/opt/llvm/lib/clang/4.0.0/include)

include_directories(${source_dir}/include) #gtest
include_directories(${PROJECT_SOURCE_DIR}/src/Core/include 
		    ${PROJECT_SOURCE_DIR}/src/Base/include)

#link_directories(/usr/local/opt/llvm/lib)
include_directories(${GAUSS_INCLUDE_DIR})
link_directories(${GAUSS_LINK_DIR})

add_subdirectory(${PROJECT_SOURCE_DIR}/src/Core)

#odd ordering is because implicit solver lives in Base currently
add_subdirectory(${PROJECT_SOURCE_DIR}/src/SolversLinear) 
add_subdirectory(${PROJECT_SOURCE_DIR}/src/Optimization)

add_subdirectory(${PROJECT_SOURCE_DIR}/src/Base)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/ParticleSystem)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/FEM)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/Collisions)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/UI)

include_directories(${GAUSS_INCLUDE_DIR})
set(Gauss_EXT_LIBS ${GAUSS_LIBS} CACHE STRING "" FORCE)
set(Gauss_EXT_INCLUDE_DIRS ${GAUSS_INCLUDE_DIR} CACHE STRING " " FORCE)
set(Gauss_EXT_LIBDIR ${GAUSS_LINK_DIR} CACHE STRING "" FORCE)

if(USE_MATLAB)
    add_subdirectory(${PROJECT_SOURCE_DIR}/src/MATLAB)
endif(USE_MATLAB)

#Init GUI Stuff
include(${PROJECT_SOURCE_DIR}/src/UI/UISetup.txt)
initUI()

#Example codes
add_subdirectory(${PROJECT_SOURCE_DIR}/src/Examples)
#Disable tests for now because I'm not using them and they only build on one machine 
#add_subdirectory(${PROJECT_SOURCE_DIR}/src/Tests)

