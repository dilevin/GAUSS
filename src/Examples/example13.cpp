#include <functional>

#include <Qt3DIncludes.h>
#include <GaussIncludes.h>
#include <ParticleSystemIncludes.h>
#include <FEMIncludes.h>
#include <Newton.h>

#include <tuple>
//Any extra things I need such as constraints
#include <ConstraintFixedPoint.h>
#include <TimeStepperEulerImplicitBFGS.h>
#include <type_traits>
using namespace Gauss;
using namespace FEM;
using namespace ParticleSystem; //For Force Spring

/* Tetrahedral finite elements */

//typedef physical entities I need
//typedef scene
typedef PhysicalSystemFEM<double, NeohookeanTet> FEMLinearTets;

typedef World<double, std::tuple<FEMLinearTets *,PhysicalSystemParticleSingle<double> *>,
std::tuple<ForceSpringFEMParticle<double> *, ForceParticlesGravity<double> *>,
std::tuple<ConstraintFixedPoint<double> *> > MyWorld;
typedef TimeStepperEulerImplicitBFGS<double, AssemblerEigenSparseMatrix<double>,
AssemblerEigenVector<double> > MyTimeStepper;

typedef Scene<MyWorld, MyTimeStepper> MyScene;


void preStepCallback(MyWorld &world) {
    // This is an example callback
}


//Test newton solver by solving for static equilibrium of a bendy bar
int main(int argc, char **argv) {
    std::cout<<"Test Neohookean FEM \n";
    
    //Setup Physics
    MyWorld world;
    
    //new code -- load tetgen files
    //Eigen::MatrixXd V(4,3);
    //Eigen::MatrixXi F(1,4);
    
    //V << 0,0,0,
      //   1,0,0,
       //  0,1,0,
        // 0,0,1;
    
    //F << 0, 1, 2, 3;
    
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    
    readTetgen(V, F, dataDir()+"/meshesTetgen/Beam/Beam.node", dataDir()+"/meshesTetgen/Beam/Beam.ele");
    
    FEMLinearTets *test = new FEMLinearTets(V,F);
    world.addSystem(test);
    
    //fixDisplacementMin(world, test);
    world.finalize(); //After this all we're ready to go (clean up the interface a bit later)
    
    auto q = mapStateEigen(world);
    q.setZero();
    
    Eigen::VectorXi indices = minVertices(test, 0);
        
    Eigen::SparseMatrix<double> P = fixedPointProjectionMatrix(indices, *test,world);
    
    MyTimeStepper stepper(0.1,P);
    
    //Display
    QGuiApplication app(argc, argv);
    
    MyScene *scene = new MyScene(&world, &stepper, preStepCallback);
    GAUSSVIEW(scene);
    
    return app.exec();
    
}
