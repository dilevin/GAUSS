#include <functional>

#include <Qt3DIncludes.h>
#include <GaussIncludes.h>
#include <FEMIncludes.h>

//Any extra things I need such as constraints
#include <ConstraintFixedPoint.h>
#include <TimeStepperEigenFitSMW.h>
#include <EigenFit.h>

using namespace Gauss;
using namespace FEM;
using namespace ParticleSystem; //For Force Spring

/* Tetrahedral finite elements */

//typedef physical entities I need

//typedef scene
typedef PhysicalSystemFEM<double, NeohookeanTet> FEMLinearTets;

typedef World<double, std::tuple<FEMLinearTets *>,
std::tuple<ForceSpringFEMParticle<double> *, ForceParticlesGravity<double> *>,
std::tuple<ConstraintFixedPoint<double> *> > MyWorld;

//typedef World<double, std::tuple<FEMLinearTets *,PhysicalSystemParticleSingle<double> *>,
//                      std::tuple<ForceSpringFEMParticle<double> *>,
//                      std::tuple<ConstraintFixedPoint<double> *> > MyWorld;
typedef TimeStepperEigenFitSMW<double, AssemblerEigenSparseMatrix<double>,
AssemblerEigenVector<double> > MyTimeStepper;

typedef Scene<MyWorld, MyTimeStepper> MyScene;


//typedef TimeStepperEigenFitSI<double, AssemblerParallel<double, AssemblerEigenSparseMatrix<double> >,
//AssemblerParallel<double, AssemblerEigenVector<double> > > MyTimeStepper;

//typedef Scene<MyWorld, MyTimeStepper> MyScene;


void preStepCallback(MyWorld &world) {
    // This is an example callback
}

int main(int argc, char **argv) {
    std::cout<<"Test Neohookean FEM \n";
    
    //Setup Physics
    MyWorld world;
    
    //new code -- load tetgen files
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    
    readTetgen(V, F, dataDir()+"/meshesTetgen/bar/barmesh5.node", dataDir()+"/meshesTetgen/bar/barmesh5.ele");

    //new code -- load tetgen files
    Eigen::MatrixXd Vf;
    Eigen::MatrixXi Ff;
    
    readTetgen(Vf, Ff, dataDir()+"/meshesTetgen/bar/barmesh8.node", dataDir()+"/meshesTetgen/bar/barmesh8.ele");

    
    
    EigenFit *test = new EigenFit(Vf,Ff,Vf,Ff);

    world.addSystem(test);
    fixDisplacementMin(world, test);
    world.finalize(); //After this all we're ready to go (clean up the interface a bit later)
    
     auto q = mapStateEigen(world);
     q.setZero();
    
    
    Eigen::VectorXi indices = minVertices(test, 0);
    Eigen::SparseMatrix<double> P = fixedPointProjectionMatrix(indices, *test,world);
    
    
     MyTimeStepper stepper(0.01,P);
    
    //Display
    QGuiApplication app(argc, argv);
    
     MyScene *scene = new MyScene(&world, &stepper, preStepCallback);
     GAUSSVIEW(scene);
    
    return app.exec();
}
