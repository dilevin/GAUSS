
cmake_minimum_required(VERSION 3.2)

project (ParticleSystem)

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${GAUSS_BINARY_DIR}/bin )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${GAUSS_BINARY_DIR}/lib )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${GAUSS_BINARY_DIR}/lib )

file(GLOB HEADER_CODE ${ParticleSystem_SOURCE_DIR}/include/*.h)
file(GLOB SRC_CODE ${ParticleSystem_SOURCE_DIR}/cpp/*.cpp)

set(GAUSS_INCLUDE_DIR ${GAUSS_INCLUDE_DIR}  ${ParticleSystem_SOURCE_DIR}/include CACHE INTERNAL "" FORCE)

include_directories(${GAUSS_INCLUDE_DIR})
add_library(ParticleSystem ${HEADER_CODE} ${SRC_CODE})
