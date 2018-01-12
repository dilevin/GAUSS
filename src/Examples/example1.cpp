#include <GaussIncludes.h>
#include <Qt3DIncludes.h>

#include <PhysicalSystemParticles.h>
#include <TimeStepperEulerImplicitLinear.h>
#include <ForceSpring.h>
#include <ConstraintFixedPoint.h>

using namespace Gauss;
using namespace ParticleSystem;

/* Particle Systems */

typedef World<double,
              std::tuple<PhysicalSystemParticleSingle<double> *>,
              std::tuple<ForceSpringParticles<double> *>,
              std::tuple<ConstraintFixedPoint<double> *> >  MyWorld;

typedef TimeStepperEulerImplicitLinear<double,
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
    
    
    ForceSpringParticles<double> *forceSpring = new ForceSpringParticles<double>(PosParticle<double>(&test->getQ()),
                                                                                 PosParticle<double>(&test1->getQ()),
                                                                                 5.0, 200.0);
    
    world.addSystem(test);
    world.addSystem(test1);
    world.addForce(forceSpring);
    world.finalize(); //After this all we're ready to go (clean up the interface a bit later)
    
    MyTimeStepper stepper(0.001);
    
    
    auto q = mapStateEigen(world);
    q.setZero();
    
    auto q1 = mapDOFEigen(test1->getQ(), world);
    q1[0] = 10.0;
    
   
    
    //Start Display
    QGuiApplication app(argc, argv);

	//Setup Scene
	MyScene *scene = new MyScene(&world, &stepper);

    GAUSSVIEW(scene);
    gaussView.startScene();
    
    return app.exec();
}

