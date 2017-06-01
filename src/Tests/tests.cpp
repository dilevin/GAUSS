#include <iostream>
/*#include <gtest/gtest.h>

//Gauss Includes
#include <PhysicalSystem.h>
#include <MultiVector.h>
#include <PhysicalSystemParticles.h>
#include <World.h>
#include <Assembler.h>
#include <Utilities.h>
#include <UtilitiesEigen.h>
#include <Assembler.h>
#include <TimeStepperEulerImplicitLinear.h>
#include <ForceSpring.h>

//FEM Stuff
#include <Element.h>
#include <PhysicalSystemFEM.h>
#include <FEMIncludes.h>

//Eigen
#include <Eigen/Dense>
#include <Eigen/Sparse>


using namespace Gauss;
using namespace ParticleSystem;


//Single Particle Tests
TEST (SingleParticleTest, Tests) {
    World<double, std::tuple<PhysicalSystemParticleSingle<double> *>, std::tuple<ForceSpring<double> *> > world;
    PhysicalSystemParticleSingle<double> *test =
    new PhysicalSystemParticleSingle<double>();
    world.addSystem(test);
    world.finalize();
    
    
    //check right amount of memory is reported and that position and velocity DOFS are in the right place
    //single particle global id of q and qdot = 0, state offset = 3
    //make sure world is taking up the correct pointer
    ASSERT_EQ(world.getSystemList().get<0>(0)->getQ().getGlobalId(), 0);
    ASSERT_EQ(world.getSystemList().get<0>(0)->getQDot().getGlobalId(), 0);
    
    //set State using offset state variable pointer
    double q[] = {0.2,5.6,4.2};
    double qDot[] = {7.6,3.1,2.9};
    
    ASSERT_EQ(3, world.getState().stateSize(0));
    ASSERT_EQ(3, world.getState().stateSize(1));
    
    world.getState().setState<0>(q, 3);
    world.getState().setState<1>(qDot, 3);
    
    
    ASSERT_EQ(world.getState()[0],0.2);
    ASSERT_EQ(world.getState()[1],5.6);
    ASSERT_EQ(world.getState()[2],4.2);
    
    ASSERT_EQ(world.getState()[3],7.6);
    ASSERT_EQ(world.getState()[4],3.1);
    ASSERT_EQ(world.getState()[5],2.9);
    
    std::pair<double *, unsigned int> ptr = world.getState().getStatePtr();

    ASSERT_EQ(6, ptr.second);
    
    Eigen::Map<Eigen::VectorXd> eigenMap(ptr.first, ptr.second);
    
    
    //Test Setting and reading state with an Eigen Map
    Eigen::VectorXd testVec(6);
    testVec[0] = 1.99;
    testVec[1] = 67.22;
    testVec[2] = 22.44;
    testVec[3] = 34.23;
    testVec[4] = 187.26;
    testVec[5] = 4001.12;
    
    eigenMap = testVec;
    
    ASSERT_EQ(world.getState()[0],1.99);
    ASSERT_EQ(world.getState()[1],67.22);
    ASSERT_EQ(world.getState()[2],22.44);
    
    ASSERT_EQ(world.getState()[3],34.23);
    ASSERT_EQ(world.getState()[4],187.26);
    ASSERT_EQ(world.getState()[5],4001.12);

    //Test Setting and reading state through dofs
    std::pair<double *, unsigned int> qPtr = world.getSystemList().get<0>(0)->getQ().getPtr(world.getState());
    std::pair<double *, unsigned int> qDotPtr = world.getSystemList().get<0>(0)->getQDot().getPtr(world.getState());
    
    //set position DOF and check
    Eigen::Map<Eigen::Vector3d> qMap(qPtr.first, qPtr.second);
    Eigen::Map<Eigen::Vector3d> qDotMap(qDotPtr.first, qDotPtr.second);
    
    qMap[0] = 4534.4;
    qMap[1] = 43573.1;
    qMap[2] = 3570.76;
    
    ASSERT_EQ(4534.4, world.getState()[0]);
    ASSERT_EQ(43573.1, world.getState()[1]);
    ASSERT_EQ(3570.76, world.getState()[2]);
    
    ASSERT_EQ(34.23, world.getState()[3]);
    ASSERT_EQ(187.26, world.getState()[4]);
    ASSERT_EQ(4001.12, world.getState()[5]);
    
    qDotMap[0] = 47384.4;
    qDotMap[1] = 423473.1;
    qDotMap[2] = 2424.76;

    //set Velocity DOF and check
    ASSERT_EQ(4534.4, world.getState()[0]);
    ASSERT_EQ(43573.1, world.getState()[1]);
    ASSERT_EQ(3570.76, world.getState()[2]);
    
    ASSERT_EQ(47384.4, world.getState()[3]);
    ASSERT_EQ(423473.1, world.getState()[4]);
    ASSERT_EQ(2424.76, world.getState()[5]);

    //Test the assembler
    Assembler<double, AssemblerImplEigenSparseMatrix> testAssembler;
    Eigen::SparseMatrix<double> testMat;
    Eigen::SparseMatrix<double> resultMat(3,3);
    test->getImpl().setMass(22.0);
    
    resultMat.setIdentity();
    resultMat *= test->getImpl().getMass();
    
    ASSEMBLEMAT(world, testAssembler, getNumQDotDOFs, getNumQDotDOFs,getMassMatrix);
    
    //testAssembler.init(world.getNumQDotDOFs(), world.getNumQDotDOFs());
    //forEach(world.getSystemList(), [&world, &testAssembler](auto a) {
      //  a->getMassMatrix(testAssembler, world.getState());
    //});
    
    testAssembler.finalize();

    ASSERT_EQ(resultMat.valuePtr()[0], testAssembler.getImpl().getMatrix().valuePtr()[0]);
    ASSERT_EQ(resultMat.valuePtr()[1], testAssembler.getImpl().getMatrix().valuePtr()[1]);
    ASSERT_EQ(resultMat.valuePtr()[2], testAssembler.getImpl().getMatrix().valuePtr()[2]);
    
    
}

TEST (SingleParticleTest, ImplicitIntegratorTest) {
    
    //This is broken for now, I pulled gravity from the integrator
    //Start at (0, 5.9,0) with initial velocity (0, -1, 0) check that after one second we're at (0,0,0)
    World<double, std::tuple<PhysicalSystemParticleSingle<double> *>, std::tuple<ForceSpring<double> *> > world;
    PhysicalSystemParticleSingle<double> *test =
    new PhysicalSystemParticleSingle<double>();
    world.addSystem(test);
    world.finalize();

    //Initial State
    Eigen::Map<Eigen::VectorXd> state = mapStateEigen(world);
    state.setZero();
    state[1] = 5.9;
    state[4] = -1.0;
    
    //Integrate
    TimeStepperEulerImplictLinear<double, Assembler<double, AssemblerImplEigenSparseMatrix>,Assembler<double, AssemblerImplEigenVector> > stepper(0.01);
    for(unsigned int ii=0; ii<10; ++ii) {
        stepper.step(world);
    }
    
    //Check
    ASSERT_EQ(0, state[0]);
    ASSERT_LE(fabs(state[1]-5.7461), 1e-5);
    ASSERT_EQ(0, state[2]);
    
    ASSERT_EQ(0, state[3]);
    ASSERT_LE(fabs(state[4]-(-1.98)), 1e-5); //This is failing currently because I removed gravity from the integrator (I want to replace it with a force)
    ASSERT_EQ(0, state[5]);
}

TEST(FEM, InitLinearFEM) {
    
    //std::array<DOFBase *, 8> dofs;
    //FEM::Element<double, 4, FEM::QuadratureExact, FEM::EnergyKineticNonLumped, FEM::EnergyPotentialNone, FEM::ShapeFunctionLinearTet> femElement(NULL, NULL, dofs);
 
    //Single tetrahedra mesh
    Eigen::MatrixXd V(4,3);
    Eigen::MatrixXi F(1,4);
 
    V << 1.00, 1.00, 1.00,
         2.00, 1.00, 1.00,
         1.00, 2.00, 1.00,
         1.00, 1.00, 2.00;
    
    F<<1, 2, 3, 0;
    
    FEM::PhysicalSystemFEM<double, LinearTet> fem(V, F);
    
    
    
}*/

int main(int argc, char **argv) {
    std::cout<<"Start Tests ..... \n";
    
   // ::testing::InitGoogleTest(&argc, argv);
    //return RUN_ALL_TESTS();
    
    //return 1;
}
