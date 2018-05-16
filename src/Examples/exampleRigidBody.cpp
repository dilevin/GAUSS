#include <functional>

#include <Qt3DIncludes.h>
#include <GaussIncludes.h>
#include <PhysicalSystemRigidBody.h>

//Any extra things I need such as constraints
#include <ConstraintFixedPoint.h>
#include <TimeStepperEulerImplicitLinear.h>

using namespace Gauss;
using namespace FEM;
using namespace RigidBodies;
using namespace ParticleSystem; //For Force Spring

/* Tetrahedral finite elements */

//typedef physical entities I need

//typedef scene
typedef PhysicalSystemRigidBody<double> RigidBody;

typedef World<double, std::tuple<RigidBody *>,
std::tuple<ForceSpringFEMParticle<double> *>,
std::tuple<ConstraintFixedPoint<double> *> > MyWorld;
typedef TimeStepperEulerImplicitLinear<double, AssemblerEigenSparseMatrix<double>,
AssemblerEigenVector<double> > MyTimeStepper;

typedef Scene<MyWorld, MyTimeStepper> MyScene;


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
    
    readTetgen(V, F, dataDir()+"/meshesTetgen/Beam/Beam.node", dataDir()+"/meshesTetgen/Beam/Beam.ele");
    
    RigidBody *test = new RigidBody(V,F);
    
    world.addSystem(test);
    //fixDisplacementMin(world, test);
    world.finalize(); //After this all we're ready to go (clean up the interface a bit later)
    
    auto q = mapStateEigen(world);
    q.setZero();
    
    //MyTimeStepper stepper(0.01);
    
    //Display
    //QGuiApplication app(argc, argv);
    
    //MyScene *scene = new MyScene(&world, &stepper, preStepCallback);
    //GAUSSVIEW(scene);
    
    //return app.exec();
}
