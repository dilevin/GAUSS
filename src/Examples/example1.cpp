#define DEBUG //Debug mode

#include <GaussIncludes.h>
#include <Qt3DIncludes.h>

#include <PhysicalSystemParticles.h>
#include <TimeStepperEulerImplicitLinear.h>
#include <ForceSpring.h>
#include <ConstraintFixedPoint.h>

using namespace Gauss;
using namespace ParticleSystem;

//1. Final do to, figure out external forces
//2. Add external spring force between DOF and a point in space
typedef World<double,
              std::tuple<PhysicalSystemParticleSingle<double> *>,
              std::tuple<ForceSpring<double> *>,
              std::tuple<ConstraintFixedPoint<double> *> >  MyWorld;

typedef TimeStepperEulerImplictLinear<double,
                                      Assembler<double, AssemblerImplEigenSparseMatrix>,
                                      Assembler<double, AssemblerImplEigenVector> > MyTimeStepper;

typedef Scene<MyWorld, MyTimeStepper> MyScene;


int main(int argc, char **argv)
{
    std::cout<<"Example 1\n";

    //Setup Physics
    MyWorld world;
    PhysicalSystemParticleSingle<double> *test =
    new PhysicalSystemParticleSingle<double>();
    
    PhysicalSystemParticleSingle<double> *test1 =
    new PhysicalSystemParticleSingle<double>();
    
    
    ForceSpring<double> *forceSpring = new ForceSpring<double>(&test->getQ(), &test1->getQ(), 5.0, 2.0);
    
    world.addSystem(test);
    world.addSystem(test1);
    world.addForce(forceSpring);
    world.finalize(); //After this all we're ready to go (clean up the interface a bit later)
    
    MyTimeStepper stepper(0.001);
    
    
    auto q = mapStateEigen(world);
    q.setZero();
    
    auto q1 = mapDOFEigen(test1->getQ(), world);
    q1[0] = 10.0;
    
    //Setup Scene
    MyScene *scene = new MyScene(&world, &stepper);
    
    //Start Display
    QGuiApplication app(argc, argv);
  
    GAUSSVIEW(scene);
    gaussView.startScene();
    
    return app.exec();
}

