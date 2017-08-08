cmake_minimum_required(VERSION 3.5)

project(Example1)

add_executable(Example1 example1.cpp)
add_executable(Example2 example2.cpp)
add_executable(Example3 example3.cpp)
if(USE_OPENMP)
	add_executable(Example4 example4.cpp)
endif(USE_OPENMP)

add_executable(Example5 example5.cpp)

target_link_libraries(Example1 Core Base UI ${GAUSS_LIBS})
target_link_libraries(Example2 Core Base FEM UI ${GAUSS_LIBS})
target_link_libraries(Example3 Core Base FEM UI ${GAUSS_LIBS})
if(USE_OPENMP)
	target_link_libraries(Example4 Core Base FEM UI ${GAUSS_LIBS})
endif(USE_OPENMP)

target_link_libraries(Example5 Core Base FEM UI ${GAUSS_LIBS})

