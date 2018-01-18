#include <functional>

#include <Qt3DIncludes.h>
#include <GaussIncludes.h>
#include <ParticleSystemIncludes.h>
#include <FEMIncludes.h>
#include <tuple>
//Any extra things I need such as constraints
#include <ConstraintFixedPoint.h>
#include <TimeStepperEulerImplicitLinear.h>
#include <TimeStepperEulerImplicit.h>
#include <type_traits>
#include <UtilitiesFEM.h>
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
typedef TimeStepperEulerImplicitLinear<double, AssemblerEigenSparseMatrix<double>,
AssemblerEigenVector<double> > MyTimeStepper;

typedef Scene<MyWorld, MyTimeStepper> MyScene;


std::vector<ConstraintFixedPoint<double> *> movingConstraints;

void preStepCallback(MyWorld &world) {
    // This is an example callback

    //script some motion
    for(unsigned int jj=0; jj<movingConstraints.size(); ++jj) {
        if(movingConstraints[jj]->getImpl().getFixedPoint()[0] > -3) {
            movingConstraints[jj]->getImpl().setFixedPoint(movingConstraints[jj]->getImpl().getFixedPoint() + Eigen::Vector3d(-0.1,0,0));
        }
    }
}

int main(int argc, char **argv) {
    std::cout<<"Test Neohookean FEM \n";

    //Setup Physics
    MyWorld world;

    //new code -- load tetgen files
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    readTetgen(V, F, dataDir()+"/meshesTetgen/Beam/Beam.node", dataDir()+"/meshesTetgen/Beam/Beam.ele");

    FEMLinearTets *test = new FEMLinearTets(V,F);

    world.addSystem(test);

    fixDisplacementMin(world, test); //fix one side
    Eigen::VectorXi movingVerts = maxVertices(test, 0);//indices for moving parts

    for(unsigned int ii=0; ii<movingVerts.rows(); ++ii) {
        movingConstraints.push_back(new ConstraintFixedPoint<double>(&test->getQ()[movingVerts[ii]], Eigen::Vector3d(0,0,0)));
        world.addConstraint(movingConstraints[ii]);
    }

    world.finalize(); //After this all we're ready to go (clean up the interface a bit later)

    auto q = mapStateEigen(world);
    q.setZero();

    MyTimeStepper stepper(0.01);

    //Display
    QGuiApplication app(argc, argv);

    MyScene *scene = new MyScene(&world, &stepper, preStepCallback);
    GAUSSVIEW(scene);

    return app.exec();

}
