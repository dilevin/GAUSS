cmake_minimum_required(VERSION 3.2)

project(Tests)

add_executable(Tests tests.cpp)

target_link_libraries(Tests Core Base libgtest ${GAUSS_LIBS})
