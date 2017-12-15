# - Try to find Gauss lib
#
# This module supports requiring a minimum version, e.g. you can do
#   find_package(Gauss 3.1.2)
# to require version 3.1.2 or newer of Gauss.
#
# Once done this will define
#
#  Gauss_FOUND - system has eigen lib with correct version
#  Gauss_INCLUDE_DIR - the eigen include directory
#  Gauss_VERSION - eigen version

# Copyright (c) 2006, 2007 Montel Laurent, <montel@kde.org>
# Copyright (c) 2008, 2009 Gael Guennebaud, <g.gael@free.fr>
# Copyright (c) 2009 Benoit Jacob <jacob.benoit.1@gmail.com>
# Redistribution and use is allowed according to the terms of the 2-clause BSD license.

macro(_Gauss_check_version)
  set(Gauss_VERSION "1.0")
endmacro(_Gauss_check_version)

#are my include directories set up ?
if (Gauss_ROOT_DIR)

message(WARNING ${GAUSS_ROOT_DIR})
  # in cache already
  _Gauss_check_version()
  set(Gauss_FOUND ${Gauss_VERSION_OK})

else (Gauss_ROOT_DIR)

  find_path(Gauss_ROOT_DIR NAMES signature_of_gauss_library
      PATHS
      ${CMAKE_INSTALL_PREFIX}/include
      ${KDE4_INCLUDE_DIR}
      PATH_SUFFIXES Gauss gauss
    )

  if(Gauss_ROOT_DIR)
    _Gauss_check_version()
  endif(Gauss_ROOT_DIR)

  
  mark_as_advanced(Gauss_ROOT_DIR)

endif(Gauss_ROOT_DIR)

message(WARNING ${Gauss_ROOT_DIR})

#setup all the variables I need to run Gauss
#Read in CMakeCache to grab important variables. Just finds cache file in directory, specifying the
#file causes a silent failure which results in no variables being imported
load_cache(${Gauss_ROOT_DIR}/build/)

#include files
set(Gauss_INCLUDE_DIRS  ${LIBIGL_INCLUDE_PATH}
                        ${EIGEN3_INCLUDE_DIR}
			            ${SolversLinear_SOURCE_DIR}/include
                        ${Optimization_SOURCE_DIR}/include 
                        ${Base_SOURCE_DIR}/include 
                        ${Core_SOURCE_DIR}/include
                        ${ParticleSystem_SOURCE_DIR}/include
                        ${FEM_SOURCE_DIR}/include
                        ${Core_SOURCE_DIR}/include
                        ${Collisions_SOURCE_DIR}/include
                        ${UI_SOURCE_DIR}/include
                        )
if(APPLE)
  if(USE_OPENMP)
          set(CMAKE_C_COMPILER ${LLVM_BIN}/clang CACHE STRING "C compiler" FORCE)
          set(CMAKE_CXX_COMPILER ${LLVM_BIN}/clang CACHE STRING "C++ compiler" FORCE)
          set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} -fopenmp)
          set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -fopenmp)
          set(CMAKE_XCODE_ATTRIBUTE_CC /usr/local/opt/llvm/bin/clang)
          set(Gauss_INCLUDE_DIRS  ${Gauss_INCLUDE_DIRS} ${LLVM_INCLUDE})
          add_definitions(-DGAUSS_OPENMP)
  endif(USE_OPENMP)
endif(APPLE)

#define the initUI macro
include(${UI_SOURCE_DIR}/UISetup.txt)

#Currently for xcode builds 
#libraries
set(Gauss_LIB_DIR_DEBUG ${Gauss_ROOT_DIR}/build/lib/Debug)
set(Gauss_LIB_DIR_RELEASE ${Gauss_ROOT_DIR}/build/lib/Release)

set(Gauss_LIBS  libBase.a
                libCore.a
                libFEM.a
                libUI.a
                libCollisions.a
                ${Gauss_EXT_LIBS})

message(WARNING "INCLUDES: " ${Gauss_INCLUDE_DIRS})
message(WARNING "DEBUG LIBS: " ${Gauss_LIBS})

