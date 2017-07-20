cmake_minimum_required(VERSION 3.5)

project(Tests)

add_executable(Tests tests.cpp)

target_link_libraries(Tests Core Base gtest_main ${GAUSS_LIBS})
