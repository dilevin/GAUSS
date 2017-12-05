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

typedef World<double, std::tuple<FEMLinearTets *,PhysicalSystemParticleSingle<double> *>,
std::tuple<ForceSpringFEMParticle<double> *, ForceParticlesGravity<double> *>,
std::tuple<ConstraintFixedPoint<double> *> > MyWorld;

typedef TimeStepperEulerImplictLinear<double, AssemblerEigenSparseMatrix<double>,
AssemblerEigenVector<double> > MyTimeStepper;

typedef Scene<MyWorld, MyTimeStepper> MyScene;

typedef CollisionDetector<double, MyWorld, CollisionFloorImpl<double> > MyCollisionDetector;

int main(int argc, char **argv) {
    std::cout<<"Collisions Test \n";
    
    MyWorld world;
    
    //new code -- load tetgen files
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    
    readTetgen(V, F, dataDir()+"/meshesTetgen/Beam/Beam.node", dataDir()+"/meshesTetgen/Beam/Beam.ele");
    
    FEMLinearTets *test = new FEMLinearTets(V,F);
    world.addSystem(test);
    world.finalize();
    
    auto q = mapStateEigen(world);
    q.setZero();
    
    MyCollisionDetector cd(world, Eigen::Vector3d(0.0,1.0,0.0), Eigen::Vector3d(0.0,0.0,0.0));
    cd.detectCollisions();
    
    //print out list of collisions
    AssemblerEigenSparseMatrix<double> collisions;
    //get mass matrix
    ASSEMBLEMATINIT(collisions, cd.getNumCollisions(), world.getNumQDotDOFs());
    cd.getGradient<AssemblerEigenSparseMatrix<double>, 0>(collisions, world.getState());
    ASSEMBLEEND(collisions);
    
}
