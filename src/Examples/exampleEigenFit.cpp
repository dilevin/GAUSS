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
typedef PhysicalSystemFEM<double, NeohookeanHFixedTet> FEMLinearTets;

typedef World<double, std::tuple<FEMLinearTets *>,
std::tuple<ForceSpringFEMParticle<double> *, ForceParticlesGravity<double> *>,
std::tuple<ConstraintFixedPoint<double> *> > MyWorld;

//typedef World<double, std::tuple<FEMLinearTets *,PhysicalSystemParticleSingle<double> *>,
//                      std::tuple<ForceSpringFEMParticle<double> *>,
//                      std::tuple<ConstraintFixedPoint<double> *> > MyWorld;
//typedef TimeStepperEigenFitSMW<double, AssemblerEigenSparseMatrix<double>, AssemblerEigenVector<double>> MyTimeStepper;
typedef TimeStepperEigenFitSMW<double, AssemblerParallel<double, AssemblerEigenSparseMatrix<double>>, AssemblerParallel<double, AssemblerEigenVector<double>> > MyTimeStepper;

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
    
//    readTetgen(V, F, dataDir()+"/meshesTetgen/bar/barmesh5.node", dataDir()+"/meshesTetgen/bar/barmesh5.ele");
//    readTetgen(V, F, dataDir()+"/meshesTetgen/Beam/Beam.node", dataDir()+"/meshesTetgen/Beam/Beam.ele");
    
    readTetgen(V, F, dataDir()+"/meshesTetgen/arma/arma_6.node", dataDir()+"/meshesTetgen/arma/arma_6.ele");
    
    //new code -- load tetgen files
    Eigen::MatrixXd Vf;
    Eigen::MatrixXi Ff;

    readTetgen(Vf, Ff, dataDir()+"/meshesTetgen/arma/arma_4.node", dataDir()+"/meshesTetgen/arma/arma_4.ele");

//    readTetgen(Vf, Ff, dataDir()+"/meshesTetgen/bar/barmesh8.node", dataDir()+"/meshesTetgen/bar/barmesh8.ele");
//    readTetgen(Vf, Ff, dataDir()+"/meshesTetgen/Beam/Beam.node", dataDir()+"/meshesTetgen/Beam/Beam.ele");

    
    
    EigenFit *test = new EigenFit(V,F,Vf,Ff);
    
    for(unsigned int iel=0; iel<test->getImpl().getF().rows(); ++iel) {
        
        test->getImpl().getElement(iel)->setParameters(1e6, 0.45);
        
    }
//    test->getImpl.getElements.setParameters(1e5,0.45);
    world.addSystem(test);
    fixDisplacementMin(world, test,2);
    world.finalize(); //After this all we're ready to go (clean up the interface a bit later)
    
     auto q = mapStateEigen(world);
     q.setZero();
    
    
    Eigen::VectorXi indices = minVertices(test, 2);
    Eigen::SparseMatrix<double> P = fixedPointProjectionMatrix(indices, *test,world);
    
    
     MyTimeStepper stepper(0.01,P);
    
    //Display
    QGuiApplication app(argc, argv);
    
     MyScene *scene = new MyScene(&world, &stepper, preStepCallback);
     GAUSSVIEW(scene);
    
    return app.exec();
}
