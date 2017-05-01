cmake_minimum_required (VERSION 3.5)
include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

project (Gauss)

set(CMAKE_CXX_STANDARD 14)

#Output directories
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib )

#Global variables to configure other projects
set(GAUSS_INCLUDE_DIR “” CACHE INTERNAL "" FORCE)
set(GAUSS_BINARY_DIR ${PROJECT_BINARY_DIR} CACHE INTERNAL "" FORCE)

#add define for data directory
add_definitions(-DGAUSS_DATA_DIR=${PROJECT_SOURCE_DIR}/data)

#grab important packages
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

#Eigen
find_package(Eigen3 REQUIRED)

#GTest (From Masha Shugrina)
# Google Test
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Download and install
ExternalProject_Add(
gtest
URL https://github.com/google/googletest/archive/release-1.7.0.zip
PREFIX ${PROJECT_SOURCE_DIR}/ThirdParty/gtest
CMAKE_ARGS ${CMAKE_ARGS}
# Disable install step
INSTALL_COMMAND ""
                )

# Create a libgtest target to be used as a dependency by test programs
add_library(libgtest IMPORTED STATIC GLOBAL)
add_dependencies(libgtest gtest)

# Set gtest properties
ExternalProject_Get_Property(gtest source_dir binary_dir)
set_target_properties(libgtest PROPERTIES
		      "IMPORTED_LOCATION" ${PROJECT_BINARY_DIR}/lib/libgtest.a
                      "IMPORTED_LINK_INTERFACE_LIBRARIES" "${CMAKE_THREAD_LIBS_INIT}"
                      )

#Individual feature directories
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories("${source_dir}/include") #gtest
include_directories(${PROJECT_SOURCE_DIR}/src/Core/include 
		    ${PROJECT_SOURCE_DIR}/src/Base/include)

add_subdirectory(${PROJECT_SOURCE_DIR}/src/Core)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/Base)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/ParticleSystem)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/FEM)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/UI)

include_directories(${GAUSS_INCLUDE_DIR})

#Init GUI Stuff
include(${PROJECT_SOURCE_DIR}/src/UI/UISetup.txt)
initUI()

#Example codes
add_subdirectory(${PROJECT_SOURCE_DIR}/src/Examples)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/tests)

