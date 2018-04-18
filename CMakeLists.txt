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
set(USE_MATLAB CACHE BOOL "Use OpenMP via homebrew llvm on OSX")

if(APPLE)
    if(USE_OPENMP)
        set(LLVM_BIN CACHE PATH "CLANG Binary Directory")
        set(LLVM_LIB CACHE STRING "CLANG Lib Directory")
        set(LLVM_INCLUDE CACHE STRING "CLANG Header Directory")
    endif(USE_OPENMP)
endif(APPLE)

if(USE_PARDISO)
    set(PARDISO_LIB CACHE STRING "Full path to Pardiso Library to use")
    if (UNIX)
        set(FORTRAN_LIB CACHE STRING "Full path to libgfortran.a")
        set(LAPACK_LIB CACHE STRING "Full path to liblapack.a")
	set(BLAS_LIB CACHE STRING "Full path to libblas.a")
    endif (UNIX)
endif(USE_PARDISO)


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

    if(WIN32) 
        set (CMAKE_C_FLAGS /openmp ${CMAKE_C_FLAGS})
        set (CMAKE_CXX_FLAGS /openmp ${CMAKE_CXX_FLAGS})
		add_definitions(-DWIN32)
    ELSE()
        set (CMAKE_C_FLAGS ${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS})
        set (CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS})
        set (CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS})
    ENDIF()
    add_definitions(-DGAUSS_OPENMP)
endif(APPLE)

if(USE_PARDISO)
    get_filename_component(PARDISO_DIR ${PARDISO_LIB} DIRECTORY)

    if(APPLE)
        set(GAUSS_LIBS ${GAUSS_LIBS} ${PARDISO_LIB} CACHE INTERNAL "") 
    elseif(UNIX)
       set(GAUSS_LIBS ${GAUSS_LIBS} ${PARDISO_LIB} ${LAPACK_LIB} ${BLAS_LIB} ${FORTRAN_LIB} quadmath CACHE INTERNAL "")
       get_filename_component(FORTRAN_DIR ${FORTRAN_LIB} DIRECTORY)
       get_filename_component(LAPACK_DIR ${LAPACK_LIB} DIRECTORY)
       get_filename_component(BLAS_DIR ${BLAS_LIB} DIRECTORY)
       set(PARDISO_DIR ${LAPACK_DIR} ${BLAS_DIR} ${FORTRAN_LIB})
    endif(APPLE)

    set(GAUSS_LINK_DIR ${GAUSS_LINK_DIR} ${PARDISO_DIR} CACHE INTERNAL "")
    add_definitions(-DGAUSS_PARDISO)
endif(USE_PARDISO)

#Output directories
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib )

#add define for data directory
add_definitions(-DGAUSS_DATA_DIR=${PROJECT_SOURCE_DIR}/data)

#grab important packages
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

#libigl include directory
set(LIBIGL_INCLUDE_PATH ${PROJECT_SOURCE_DIR}/ThirdParty/libigl/include CACHE FILEPATH "Root include directory for libigl")

#Spectra
set(SPECTRA_INCLUDE_PATH ${PROJECT_SOURCE_DIR}/ThirdParty/Spectra/include CACHE FILEPATH "Root include directory for SPECTRA")
set(GAUSS_INCLUDE_DIR ${GAUSS_INCLUDE_DIR} ${SPECTRA_INCLUDE_DIR}/include CACHE INTERNAL "")

#directories for FCL and dependencies
set(FCL_DIR ${PROJECT_SOURCE_DIR}/ThirdParty/fcl)
set(LIBCCD_DIR ${PROJECT_SOURCE_DIR}/ThirdParty/libccd)
set(OCTOMAP_DIR ${PROJECT_SOURCE_DIR}/ThirdParty/octomap)

#augment GAUSS_INCLUDE_DIRS
set(${GAUSS_INCLUDE_DIR} ${GAUSS_INCLUDE_DIR} ${FCL_DIR}/include)
set(${GAUSS_INCLUDE_DIR} ${GAUSS_INCLUDE_DIR} ${LIBCCD_DIR}/src)
set(${GAUSS_INCLUDE_DIR} ${GAUSS_INCLUDE_DIR} ${OCTOMAP_DIR}/octomap/include)

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

#Set parameters for third party libraries
#libccd
set(BUILD_DOCUMENTATION OFF CACHE INTERNAL "" FORCE)
set(BUILD_OCTOVIS_SUBPROJECT  OFF CACHE INTERNAL "" FORCE)
set(BUILD_GMOCK OFF CACHE INTERNAL "" FORCE)
set(BUILD_GTEST OFF CACHE INTERNAL ""  FORCE)
set(BUILD_DYNAMICETD3D_SUBPROJECT OFF CACHE INTERNAL "" FORCE)
set(BUILD_TESTING OFF CACHE INTERNAL FORCE "")
set(OCTOMAP_OMP OFF CACHE INTERNAL ""  FORCE)
set(FCL_BUILD_TESTS OFF CACHE INTERNAL "" FORCE)
set(FCL_WITH_OCTOMAP OFF CACHE INTERNAL "" FORCE)
set(BUILD_SHARED_LIBS OFF CACHE INTERNAL "" FORCE) 
set(FCL_STATIC_LIBRARY ON CACHE INTERNAL "" FORCE)
set(ENABLE_DOUBLE_PRECISION ON CACHE INTERNAL "" FORCE)
#set ccd include directory 
#set(octomap_DIR ${PROJECT_SOURCE_DIR}/ThirdParty/octomap/octomap CACHE INTERNAL  "" FORCE)

set(CMAKE_DISABLE_FIND_PACKAGE_ccd ON CACHE INTERNAL "" FORCE)
set(ccd_FOUND ON  CACHE INTERNAL "" FORCE) 
set(CCD_INCLUDE_DIR ${LIBCCD_DIR}/src CACHE INTERNAL "" FORCE)
set(CCD_INCLUDE_DIRS ${LIBCCD_DIR}/src CACHE INTERNAL "" FORCE)
set(CCD_LIBRARIES "${CCD_LIBRARY}" "${M_LIBRARY}")
set(ccd_FOUND ON)

#Third Party libraries
#fix for odd build behavior of libccd on windows when using static linking.
add_definitions("-DCCD_STATIC_DEFINE")
add_subdirectory(${LIBCCD_DIR})
#add_subdirectory(${OCTOMAP_DIR})
add_subdirectory(${FCL_DIR})

message(WARNING ${BUILD_GTEST})
#Main Gauss libraries
add_subdirectory(${PROJECT_SOURCE_DIR}/src/Core)

#odd ordering is because implicit solver lives in Base currently
add_subdirectory(${PROJECT_SOURCE_DIR}/src/SolversLinear) 
add_subdirectory(${PROJECT_SOURCE_DIR}/src/Optimization)

add_subdirectory(${PROJECT_SOURCE_DIR}/src/Base)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/ParticleSystem)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/FEM)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/Embeddings)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/Collisions)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/UI)

set(GAUSS_INCLUDE_DIR ${GAUSS_INCLUDE_DIR} ${PROJECT_SOURCE_DIR}/ThirdParty/LBFGS++/include)

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

