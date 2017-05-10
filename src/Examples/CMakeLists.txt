cmake_minimum_required(VERSION 3.5)

project(Example1)

add_executable(Example1 example1.cpp)
add_executable(Example2 example2.cpp)
add_executable(Example3 example3.cpp)

target_link_libraries(Example1 Core Base UI)
target_link_libraries(Example2 Core Base FEM UI)
target_link_libraries(Example3 Core Base FEM UI)
