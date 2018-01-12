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
add_executable(Example10 example10.cpp)

add_executable(Example12 example12.cpp)
target_link_libraries(Example12 Core Base Collisions FEM UI ${GAUSS_LIBS})

if(USE_PARDISO)
	add_executable(Example11 example11.cpp)
	target_link_libraries(Example11 Core Base FEM Optimization UI ${GAUSS_LIBS})
endif(USE_PARDISO)

add_executable(Example13 example13.cpp)

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
target_link_libraries(Example10 Core Base FEM UI ${GAUSS_LIBS})
target_link_libraries(Example13 Core Base FEM UI ${GAUSS_LIBS})


