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
add_executable(ExampleNeohookeanTet exampleNeohookeanTet.cpp)
add_executable(ExampleNeohookeanHFixedTet exampleNeohookeanHFixedTet.cpp)
add_executable(ExampleEigenFit exampleEigenFit.cpp)
add_executable(Example10 example10.cpp)

add_executable(ExampleCollisionsFloor exampleCollisionsFloor.cpp)
target_link_libraries(ExampleCollisionsFloor Core Base Collisions FEM UI ${GAUSS_LIBS})


add_executable(ExampleImplicitIntegration exampleImplicitIntegration.cpp)
target_link_libraries(ExampleImplicitIntegration Core Base FEM Optimization UI ${GAUSS_LIBS})


add_executable(Example13 example13.cpp)
add_executable(Example14 example14.cpp)
add_executable(Example15 example15.cpp)
add_executable(Example16 example16.cpp)
add_executable(Example17 example17.cpp)
add_executable(Example18 example18.cpp)
add_executable(ExampleTestFCL exampleTestFCL.cpp)
add_executable(ExampleTwistingBar exampleTwistingBar.cpp)
add_executable(ExampleRigidBody exampleRigidBody.cpp)
add_executable(ExampleNewmarkNonlinear exampleNewmarkNonlinear.cpp)
add_executable(ExamplePrincipalStretch examplePrincipalStretch.cpp)
add_executable(ExampleLoubignac exampleLoubignac.cpp)

#add_executable(ExampleFloatingFrame exampleFloatingFrame.cpp)

target_link_libraries(Example1 Core Base UI ${GAUSS_LIBS})
target_link_libraries(Example2 Core Base FEM UI ${GAUSS_LIBS})
target_link_libraries(Example3 Core Base FEM UI ${GAUSS_LIBS})

if(USE_OPENMP)
	target_link_libraries(Example4 Core Base FEM UI ${GAUSS_LIBS})
endif(USE_OPENMP)

add_executable(ExampleLinearModes exampleLinearModes.cpp)
target_link_libraries(ExampleLinearModes Core Base FEM UI ${GAUSS_LIBS})

target_link_libraries(Example5 Core Base FEM UI ${GAUSS_LIBS})
target_link_libraries(Example7 Core Base FEM UI ${GAUSS_LIBS})
target_link_libraries(ExampleNeohookeanTet Core Base FEM UI ${GAUSS_LIBS})
target_link_libraries(ExampleNeohookeanHFixedTet Core Base FEM UI ${GAUSS_LIBS})
target_link_libraries(ExampleEigenFit Core Base FEM UI ${GAUSS_LIBS})
target_link_libraries(Example10 Core Base FEM UI ${GAUSS_LIBS})
target_link_libraries(Example13 Core Base FEM UI ${GAUSS_LIBS})
target_link_libraries(Example14 Core Base FEM UI ${GAUSS_LIBS})
target_link_libraries(Example15 Core Base FEM UI ${GAUSS_LIBS})
target_link_libraries(Example16 Core Base FEM UI ${GAUSS_LIBS})
target_link_libraries(Example17 Core Base FEM UI ${GAUSS_LIBS})
target_link_libraries(Example18 Core Base FEM UI ${GAUSS_LIBS})
target_link_libraries(ExampleTestFCL Core Base FEM fcl ccd UI ${GAUSS_LIBS})
target_link_libraries(ExampleTwistingBar Core Base FEM UI ${GAUSS_LIBS})
target_link_libraries(ExampleRigidBody Core Base RigidBodies UI ${GAUSS_LIBS})
target_link_libraries(ExampleNewmarkNonlinear Core Base FEM UI ${GAUSS_LIBS})
target_link_libraries(ExamplePrincipalStretch Core Base FEM UI ${GAUSS_LIBS})
target_link_libraries(ExampleLoubignac igl::core igl::opengl_glfw igl::opengl_glfw_imgui Core Base FEM UI ${GAUSS_LIBS})

#target_link_libraries(ExampleFloatingFrame Core Base RigidBodies UI ${GAUSS_LIBS})
