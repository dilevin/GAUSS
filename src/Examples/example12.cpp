#include <functional>
#include <Qt3DIncludes.h>
#include <GaussIncludes.h>
#include <ParticleSystemIncludes.h>
#include <FEMIncludes.h>
#include <CollisionDetector.h>
#include <CollisionsFloor.h>

#include <tuple>
//Any extra things I need such as constraints
#include <ConstraintFixedPoint.h>
#include <TimeStepperEulerImplicitLinear.h>
#include <TimeStepperEulerImplicit.h>
#include <type_traits>
using namespace Gauss;
using namespace Collisions;
using namespace FEM;
using namespace ParticleSystem;

typedef PhysicalSystemFEM<double, NeohookeanTet> FEMLinearTets;
typedef ConstraintCollisionDetector<double, CollisionFloorImpl> MyCollisionDetector;

typedef World<double, std::tuple<FEMLinearTets *,PhysicalSystemParticleSingle<double> *>,
std::tuple<ForceSpringFEMParticle<double> *, ForceParticlesGravity<double> *>,
std::tuple<ConstraintFixedPoint<double> *, MyCollisionDetector *> > MyWorld;

typedef TimeStepperEulerImplictLinear<double, AssemblerEigenSparseMatrix<double>,
AssemblerEigenVector<double> > MyTimeStepper;

typedef Scene<MyWorld, MyTimeStepper> MyScene;

int main(int argc, char **argv) {
    std::cout<<"Collisions Test \n";
    
    MyWorld world;
    
    //new code -- load tetgen files
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    
    readTetgen(V, F, dataDir()+"/meshesTetgen/Beam/Beam.node", dataDir()+"/meshesTetgen/Beam/Beam.ele");
    
    FEMLinearTets *test = new FEMLinearTets(V,F);
    MyCollisionDetector cd(std::ref(world), Eigen::Vector3d(0.0,1.0,0.0), Eigen::Vector3d(0.0,0.0,0.0));
    world.addSystem(test);
    world.addInequalityConstraint(&cd);
    world.finalize();
    
    auto q = mapStateEigen(world);
    q.setZero();
    
    
    
    cd.update(world);
    world.updateInequalityConstraints();
    
    //print out list of collisions
    AssemblerEigenSparseMatrix<double> collisions;
    //get mass matrix
    ASSEMBLEMATINIT(collisions, cd.getNumRows(), world.getNumQDotDOFs());
    cd.getGradient<MyWorld,AssemblerEigenSparseMatrix<double>, 0>(collisions, world, world.getState());
    ASSEMBLEEND(collisions);
    
    toMatlab(*collisions, "./constraintJacobian.txt");
    
}
