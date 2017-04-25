cmake_minimum_required(VERSION 3.5)
project (UI)

include(UISetup.txt)
initUI()

file(GLOB HEADER_CODE ${UI_SOURCE_DIR}/include/*.h)
file(GLOB SRC_CODE ${UI_SOURCE_DIR}/cpp/*.cpp)

set(GAUSS_INCLUDE_DIR ${GAUSS_INCLUDE_DIR}  ${UI_SOURCE_DIR}/include CACHE INTERNAL "" FORCE)

include_directories(${GAUSS_INCLUDE_DIR})
add_library(UI ${HEADER_CODE} ${SRC_CODE})
target_link_libraries(UI ${UI_LIBS})

