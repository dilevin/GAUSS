cmake_minimum_required(VERSION 3.2)

project(Example1)

add_executable(Example1 example1.cpp)
add_executable(Example2 example2.cpp)
add_executable(Example3 example3.cpp)
if(USE_OPENMP)
	add_executable(Example4 example4.cpp)
endif(USE_OPENMP)

add_executable(Example5 example5.cpp)
add_executable(Example7 example7.cpp)
add_executable(Example9 example9.cpp)

target_link_libraries(Example1 Core Base UI ${GAUSS_LIBS})
target_link_libraries(Example2 Core Base FEM UI ${GAUSS_LIBS})
target_link_libraries(Example3 Core Base FEM UI ${GAUSS_LIBS})
if(USE_OPENMP)
	target_link_libraries(Example4 Core Base FEM UI ${GAUSS_LIBS})
endif(USE_OPENMP)

if(USE_SPECTRA)
	add_executable(Example6 example6.cpp)
	target_link_libraries(Example6 Core Base FEM UI ${GAUSS_LIBS})
endif(USE_SPECTRA)

target_link_libraries(Example5 Core Base FEM UI ${GAUSS_LIBS})
target_link_libraries(Example7 Core Base FEM UI ${GAUSS_LIBS})
target_link_libraries(Example9 Core Base FEM UI ${GAUSS_LIBS})


