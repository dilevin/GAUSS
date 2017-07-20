cmake_minimum_required (VERSION 3.5)
include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

project (Gauss)

set(CMAKE_CXX_STANDARD 14)

#this code uses OpenMP which is supported almost everywhere except for OS X clang
#to get around this you can choose to install clang via homebrew to get access to OpenMP

#expose USE_OPENMP variable for OS X only
if(APPLE)
    if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        #message(WARNING "Here")
        set(USE_OPENMP CACHE BOOL "Use OpenMP via homebrew llvm on OSX")
    endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
endif(APPLE)

#do you want to use Pardiso
set(USE_PARDISO CACHE BOOL "Use Pardiso if available on system")

if(APPLE)
    if(USE_OPENMP)
        set(LLVM_BIN CACHE STRING "CLANG Binary Directory")
        set(LLVM_LIB CACHE STRING "CLANG Lib Directory")
        set(LLVM_INCLUDE CACHE STRING "CLANG Header Directory")
    endif(USE_OPENMP)
endif(APPLE)

if(USE_PARDISO)
    set(PARDISO_LIB CACHE STRING "Pardiso Library to use")
endif(USE_PARDISO)


#find_package(OpenMP REQUIRED)
#Global variables to configure other projects
set(GAUSS_INCLUDE_DIR  CACHE INTERNAL "")
set(GAUSS_LINK_DIR  CACHE INTERNAL "")
set(GAUSS_LIBS  CACHE INTERNAL "")
set(GAUSS_BINARY_DIR ${PROJECT_BINARY_DIR} CACHE INTERNAL "")

if(APPLE)
    if(USE_OPENMP)
        set(CMAKE_C_COMPILER ${LLVM_BIN}/clang CACHE STRING "C compiler" FORCE)
        set(CMAKE_CXX_COMPILER ${LLVM_BIN}/clang CACHE STRING "C++ compiler" FORCE)
        set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} -fopenmp)
        set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -fopenmp)
        set(CMAKE_EXE_LINKER_FLAGS -lc++)
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

#Output directories
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib )

#add define for data directory
add_definitions(-DGAUSS_DATA_DIR=${PROJECT_SOURCE_DIR}/data)

#grab important packages
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

#libigl include directory
set(LIBIGL_INCLUDE_PATH "path" CACHE FILEPATH "Root include directory for libigl")

#Eigen
find_package(Eigen3 REQUIRED)

#GTest (From Masha Shugrina)
### Google Test setup
# Download and unpack googletest at configure time
configure_file(CMakeLists.txt.gtest googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_BINARY_DIR}/googletest-src
                 ${CMAKE_BINARY_DIR}/googletest-build)

# The gtest/gtest_main targets carry header search path
# dependencies automatically when using CMake 2.8.11 or
# later. Otherwise we have to add them here ourselves.
if (CMAKE_VERSION VERSION_LESS 2.8.11)
  include_directories("${gtest_SOURCE_DIR}/include")
endif()
### End of Google Test setup


#Individual feature directories
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${LIBIGL_INCLUDE_PATH})
#include_directories(/usr/local/opt/llvm/lib/clang/4.0.0/include)

include_directories("${source_dir}/include") #gtest
include_directories(${PROJECT_SOURCE_DIR}/src/Core/include 
		    ${PROJECT_SOURCE_DIR}/src/Base/include)

#link_directories(/usr/local/opt/llvm/lib)
include_directories(${GAUSS_INCLUDE_DIR})
link_directories(${GAUSS_LINK_DIR})

add_subdirectory(${PROJECT_SOURCE_DIR}/src/Core)

#odd ordering is because implicit solver lives in Base currently
add_subdirectory(${PROJECT_SOURCE_DIR}/src/SolversLinear) 


add_subdirectory(${PROJECT_SOURCE_DIR}/src/Base)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/ParticleSystem)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/FEM)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/UI)

include_directories(${GAUSS_INCLUDE_DIR})
set(Gauss_EXT_LIBS ${GAUSS_LIBS} CACHE STRING "")

#Init GUI Stuff
include(${PROJECT_SOURCE_DIR}/src/UI/UISetup.txt)
initUI()

#Example codes
add_subdirectory(${PROJECT_SOURCE_DIR}/src/Examples)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/tests)

