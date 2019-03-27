#include <functional>

#include <Qt3DIncludes.h>
#include <GaussIncludes.h>
#include <FEMIncludes.h>
#include <igl/readDMAT.h>

//Any extra things I need such as constraints
#include <ConstraintFixedPoint.h>
#include <TimeStepperEulerImplicit.h>
#include <TimeStepperStatic.h>

using namespace Gauss;
using namespace FEM;
using namespace ParticleSystem; //For Force Spring

/* Tetrahedral finite elements */

//typedef physical entities I need

//typedef scene
typedef PhysicalSystemFEM<double, MuscleTet> FEMLinearTets;

typedef World<double, std::tuple<FEMLinearTets *,PhysicalSystemParticleSingle<double> *>,
                      std::tuple<ForceSpringFEMParticle<double> *>,
                      std::tuple<ConstraintFixedPoint<double> *> > MyWorld;
typedef TimeStepperEulerImplicit<double, AssemblerEigenSparseMatrix<double>,
AssemblerEigenVector<double> > MyTimeStepper;

typedef Scene<MyWorld, MyTimeStepper> MyScene;
FEMLinearTets *test;
double muscleStart = 10000;
Eigen::MatrixXd Ufibre;
Eigen::VectorXi imuscle;

void preStepCallback(MyWorld &world) {
    // This is an example callback
    // if(muscleStart > 1.07e8){
        // muscleStart = 100;
    // }
    // muscleStart = 1.05*muscleStart;
    std::cout<<muscleStart<<std::endl;
    // for(int m=0; m<imuscle.size(); m++){
    //     Eigen::Vector3d uvec = Ufibre.row(imuscle[m]);
    //     test->getImpl().getElements()[imuscle[m]]->setMuscleParameters(muscleStart, uvec);
    // }
}

int main(int argc, char **argv) {
    std::cout<<"Test Neohookean FEM \n";
    
    //Setup Physics
    MyWorld world;
    
    //new code -- load tetgen files
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    
    igl::readDMAT("/home/vismay/recode/fast_muscles/data/simple_muscle/generated_files/tet_mesh_T.dmat", F);
    igl::readDMAT("/home/vismay/recode/fast_muscles/data/simple_muscle/generated_files/tet_mesh_V.dmat", V);
    igl::readDMAT("/home/vismay/recode/fast_muscles/data/simple_muscle/generated_files/combined_fiber_directions.dmat", Ufibre);
    igl::readDMAT("/home/vismay/recode/fast_muscles/data/simple_muscle/generated_files/muscle_muscle_indices.dmat", imuscle);
    // readTetgen(V, F, dataDir()+"/meshesTetgen/Beam/Beam.node", dataDir()+"/meshesTetgen/Beam/Beam.ele");

    test = new FEMLinearTets(V,F);

    Eigen::Vector3d x(0,1,0);
    Eigen::Vector3d zer(0,0,0);
    for(auto element: test->getImpl().getElements()) {
        element->EnergyKineticNonLumped::setDensity(100.0);//1000.0);
        element->setParameters(1e5, 0.49);
        element->setMuscleParameters(0,x);
        element->setGravity(zer);
    }
    for(int m=0; m<imuscle.size(); m++){
        Eigen::Vector3d uvec = Ufibre.row(imuscle[m]);
        test->getImpl().getElements()[imuscle[m]]->EnergyKineticNonLumped::setDensity(10000);
        test->getImpl().getElements()[imuscle[m]]->setParameters(60000, 0.49);
        test->getImpl().getElements()[imuscle[m]]->setMuscleParameters(muscleStart, uvec);
    }




    world.addSystem(test);
    fixDisplacementMax(world, test, 1);
    world.finalize(); //After this all we're ready to go (clean up the interface a bit later)
    
    auto q = mapStateEigen(world);
    q.setZero();
    
    MyTimeStepper stepper(0.1,10);
    
    //Display
    QGuiApplication app(argc, argv);
    
    MyScene *scene = new MyScene(&world, &stepper, preStepCallback);
    GAUSSVIEW(scene);
    
    return app.exec();
}