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
#include <igl/slice.h>
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
Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::VectorXi movingVerts;

int current_frame = 0;

void preStepCallback(MyWorld &world) {
    current_frame++;
    // This is an example callback

    //script some motion
    Eigen::AngleAxis<double> rot(0.01 * current_frame, Eigen::Vector3d(1.0,0.0,0.0));

    for(unsigned int jj=0; jj<movingConstraints.size(); ++jj) {
        if(movingConstraints[jj]->getImpl().getFixedPoint()[0] > -3) {
            Eigen::Vector3d v = V.row(movingVerts[jj]);
            Eigen::Vector3d new_u = rot * v - v - mapDOFEigen(movingConstraints[jj]->getDOF(0), world.getState());
            movingConstraints[jj]->getImpl().setFixedPoint(rot*v, new_u/0.01);
        }
    }
}

int main(int argc, char **argv) {
    std::cout<<"Test Neohookean FEM \n";

    //Setup Physics
    MyWorld world;

    //new code -- load tetgen files

    Gauss::Optimization::doSomething([](){Gauss::Optimization::printVal(1.5);});
    /*readTetgen(V, F, dataDir()+"/meshesTetgen/Beam/Beam.node", dataDir()+"/meshesTetgen/Beam/Beam.ele");

    FEMLinearTets *test = new FEMLinearTets(V,F);

    world.addSystem(test);

    fixDisplacementMin(world, test); //fix one side
    movingVerts = maxVertices(test, 0);//indices for moving parts

    for(unsigned int ii=0; ii<movingVerts.rows(); ++ii) {
        movingConstraints.push_back(new ConstraintFixedPoint<double>(&test->getQ()[movingVerts[ii]], Eigen::Vector3d(0,0,0)));
        world.addConstraint(movingConstraints[ii]);
    }

    world.finalize(); //After this all we're ready to go (clean up the interface a bit later)

    auto q = mapStateEigen(world);
    q.setZero();

    MyTimeStepper stepper(0.01, 1);

    //Display
    QGuiApplication app(argc, argv);

    MyScene *scene = new MyScene(&world, &stepper, preStepCallback);
    GAUSSVIEW(scene);

    return app.exec();*/

}
