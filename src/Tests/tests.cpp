#include <iostream>
#include <gtest/gtest.h>

//Gauss Includes
#include <PhysicalSystem.h>
#include <MultiVector.h>
#include <PhysicalSystemParticles.h>
#include <World.h>
#include <Assembler.h>
#include <AssemblerParallel.h>
#include <Utilities.h>
#include <UtilitiesEigen.h>
#include <Assembler.h>
#include <AssemblerMVP.h>
#include <TimeStepperEulerImplicitLinear.h>
#include <ForceSpring.h>
#include <ConstraintFixedPoint.h>

//FEM Stuff
#include <Element.h>
#include <PhysicalSystemFEM.h>
#include <FEMIncludes.h>

//Eigen
#include <Eigen/Dense>
#include <Eigen/Sparse>

//CG Solver
#include <SolverCG.h>

using namespace Gauss;
using namespace ParticleSystem;


//Single Particle Tests
/*TEST (SingleParticleTest, Tests) {
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

TEST(MVP, TestMVP) {
    
    using namespace Gauss;
    using namespace FEM;
    
    typedef PhysicalSystemFEM<double, LinearTet> FEMLinearTets;
    
    typedef World<double, std::tuple<FEMLinearTets *>, std::tuple<ForceSpring<double> *>, std::tuple<ConstraintFixedPoint<double> *> > MyWorld;
    
    
    //init FEM system
    //Setup Physics
    MyWorld world;
    
    //new code -- load tetgen files
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    
    readTetgen(V, F, dataDir()+"/meshesTetgen/Beam/Beam.node", dataDir()+"/meshesTetgen/Beam/Beam.ele");
    
    FEMLinearTets *test = new FEMLinearTets(V,F);
    
    world.addSystem(test);
    fixDisplacementMin(world, test);
    world.finalize(); //After this all we're ready to go (clean up the interface a bit later)
    
    auto q = mapStateEigen(world);
    q.setZero();

    //test using stiffness matrix
    //explicit assembly
    AssemblerEigenSparseMatrix<double> assembler;
    ASSEMBLEMATINIT(assembler, world.getNumQDotDOFs(), world.getNumQDotDOFs());
    ASSEMBLELIST(assembler, world.getSystemList(), getStiffnessMatrix);
    ASSEMBLELIST(assembler, world.getForceList(), getStiffnessMatrix);
    ASSEMBLEEND(assembler);
    
    AssemblerMVPEigen<double> assemblerMVP;
    Eigen::VectorXd bAssembled;
    Eigen::VectorXd x;
    
    for(unsigned int ii=0; ii<100; ++ii) {
        x = 10.0*Eigen::MatrixXd::Random(world.getNumQDotDOFs(),1);
    
        Eigen::VectorXd bAssembled = (*assembler)*x;
        
        //compare mvp evaluated with random vector to assembly free MVP
        ASSEMBLEMATINIT(assemblerMVP, world.getNumQDotDOFs(), world.getNumQDotDOFs());
        assemblerMVP.getImpl().setX(x);
        
        //multiplying is happening during assembly (neat huh ? )
        ASSEMBLELIST(assemblerMVP, world.getSystemList(), getStiffnessMatrix);
        ASSEMBLELIST(assemblerMVP, world.getForceList(), getStiffnessMatrix);
        ASSEMBLEEND(assemblerMVP);
        
        double tol = 1e-8; //we'll consider this close enough
        
        ASSERT_LE((bAssembled - *assemblerMVP).norm() / x.norm(), tol);
    }
    
    
}

// TODO(lawson 19/07/17): Implement version without PARDISO?
#ifdef GAUSS_PARDISO
TEST(MVP, TestCG) {
    
    //Test CG vs. Direct Solve
    
    //Form implicit integration matrix solve against some random vectors, check CG vs direct solve
    using namespace FEM;
    
    typedef PhysicalSystemFEM<double, LinearTet> FEMLinearTets;
    
    typedef World<double, std::tuple<FEMLinearTets *>, std::tuple<ForceSpring<double> *>, std::tuple<ConstraintFixedPoint<double> *> > MyWorld;
    
    
    //init FEM system
    //Setup Physics
    MyWorld world;
    
    //new code -- load tetgen files
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    
    readTetgen(V, F, dataDir()+"/meshesTetgen/Beam/Beam.node", dataDir()+"/meshesTetgen/Beam/Beam.ele");
    
    FEMLinearTets *test = new FEMLinearTets(V,F);
    
    world.addSystem(test);
    fixDisplacementMin(world, test);
    world.finalize(); //After this all we're ready to go (clean up the interface a bit later)
    
    auto q = mapStateEigen(world);
    q.setZero();
    
    //test using stiffness matrix
    //explicit assembly
    AssemblerEigenSparseMatrix<double> M, K;
    ASSEMBLEMATINIT(M, world.getNumQDotDOFs(), world.getNumQDotDOFs());
    ASSEMBLELIST(M, world.getSystemList(), getMassMatrix);
    ASSEMBLEEND(M);
    
    ASSEMBLEMATINIT(K, world.getNumQDotDOFs(), world.getNumQDotDOFs());
    ASSEMBLELIST(K, world.getSystemList(), getStiffnessMatrix);
    ASSEMBLELIST(K, world.getForceList(), getStiffnessMatrix);
    ASSEMBLEEND(K);
    
    //random force vector
    double dt = 0.01;
    Eigen::VectorXd f;
    Eigen::VectorXd x;
    Eigen::SparseMatrix<double, Eigen::RowMajor> A = ((*M) + dt*dt*(*K)).eval();
    f = 10.0*Eigen::MatrixXd::Random(world.getNumQDotDOFs(),1);
    x = 0*f;
    
    SolverPardiso<Eigen::SparseMatrix<double, Eigen::RowMajor> > direct;
    direct.solve(A,f);
    
    SolverCG<double, Eigen::VectorXd> pcg(1e-8);
    
    std::cout<<"CG\n";
    pcg.solve(x, [&A](auto &y)->auto {return A*y;}, f, 1000000000);
    
    double tol=1e-10;
    
    ASSERT_LE((x-direct.getX()).norm()/f.norm(), tol);
}
#endif

// TODO(lawson 19/07/17): Implement version without PARDISO?
#ifdef GAUSS_PARDISO
TEST(MVP, TestCG2) {
    
    //Test Assembly free CG vs. Direct Solve
    using namespace FEM;
    
    typedef PhysicalSystemFEM<double, LinearTet> FEMLinearTets;
    
    typedef World<double, std::tuple<FEMLinearTets *>, std::tuple<ForceSpring<double> *>, std::tuple<ConstraintFixedPoint<double> *> > MyWorld;
    
    
    //init FEM system
    //Setup Physics
    MyWorld world;
    
    //new code -- load tetgen files
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    
    readTetgen(V, F, dataDir()+"/meshesTetgen/Beam/Beam.node", dataDir()+"/meshesTetgen/Beam/Beam.ele");
    
    FEMLinearTets *test = new FEMLinearTets(V,F);
    
    world.addSystem(test);
    fixDisplacementMin(world, test);
    world.finalize(); //After this all we're ready to go (clean up the interface a bit later)
    
    auto q = mapStateEigen(world);
    q.setZero();
    
    //test using stiffness matrix
    //explicit assembly
    AssemblerEigenSparseMatrix<double> M, K;
    ASSEMBLEMATINIT(M, world.getNumQDotDOFs(), world.getNumQDotDOFs());
    ASSEMBLELIST(M, world.getSystemList(), getMassMatrix);
    ASSEMBLEEND(M);
    
    ASSEMBLEMATINIT(K, world.getNumQDotDOFs(), world.getNumQDotDOFs());
    ASSEMBLELIST(K, world.getSystemList(), getStiffnessMatrix);
    ASSEMBLELIST(K, world.getForceList(), getStiffnessMatrix);
    ASSEMBLEEND(K);
    
    //random force vector
    double dt = 0.01;
    Eigen::VectorXd f;
    Eigen::VectorXd x;
    Eigen::SparseMatrix<double, Eigen::RowMajor> A = ((*M) + dt*dt*(*K)).eval();
    f = 10.0*Eigen::MatrixXd::Random(world.getNumQDotDOFs(),1);
    x = 0*f;
    
    SolverPardiso<Eigen::SparseMatrix<double, Eigen::RowMajor> > direct;
    direct.solve(A,f);
    
    SolverCG<double, Eigen::VectorXd> pcg(1e-8);
    
    
    AssemblerParallel<double, AssemblerMVPEigen<double> > Mv, Kv;
    //AssemblerMVPEigen<double> Mv, Kv;
    
    //assembly-free mvp
    auto mvp = [&world, &Mv, &Kv, &dt](auto &y) {
        ASSEMBLEMATINIT(Kv, world.getNumQDotDOFs(), world.getNumQDotDOFs());
        Kv.getImpl()*=y;
        
        //multiplying is happening during assembly (neat huh ? )
        ASSEMBLELIST(Kv, world.getSystemList(), getStiffnessMatrix);
        ASSEMBLELIST(Kv, world.getForceList(), getStiffnessMatrix);
        ASSEMBLEEND(Kv);

        ASSEMBLEMATINIT(Mv, world.getNumQDotDOFs(), world.getNumQDotDOFs());
        Mv.getImpl()*=y;
        
        //multiplying is happening during assembly (neat huh ? )
        ASSEMBLELIST(Mv, world.getSystemList(), getMassMatrix);
        ASSEMBLELIST(Mv, world.getForceList(), getMassMatrix);
        ASSEMBLEEND(Mv);
        
        return (*Mv) + dt*dt*(*Kv);
    };
    
    std::cout<<"CG\n";
    pcg.solve(x, mvp, f, 100000000);
    
    double tol=1e-10;
    
    std::cout<<f.norm()<<"\n";
    ASSERT_LE((x-direct.getX()).norm()/f.norm(), tol);

    
}
#endif

int main(int argc, char **argv) {
    std::cout<<"Start Tests ..... \n";
    
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
    
    //return 1;
}
