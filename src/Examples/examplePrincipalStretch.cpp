//Example using principal stretch-based materials
#include <functional>

#include <Qt3DIncludes.h>
#include <GaussIncludes.h>
#include <ParticleSystemIncludes.h>
#include <FEMIncludes.h>
#include <Newton.h>

#include <iostream>
#include <tuple>
//Any extra things I need such as constraints
#include <ConstraintFixedPoint.h>
#include <TimeStepperEulerImplicit.h>
#include <TimeStepperNewmark.h>
#include <TimeStepperEulerImplicitLinear.h>
#include <TimeStepperEulerImplicit.h>
#include <type_traits>

using namespace Gauss;
using namespace FEM;
using namespace ParticleSystem; //For Force Spring

//build specific principal stretch material
template<typename DataType, typename ShapeFunction>
using  EnergyPSNH = EnergyPrincipalStretch<DataType, ShapeFunction, PSNeohookean>;

template<typename DataType, typename ShapeFunction>
using  EnergyPSARAP = EnergyPrincipalStretch<DataType, ShapeFunction, PSARAP>;

template<typename DataType, typename ShapeFunction>
using  EnergyPSCoRot = EnergyPrincipalStretch<DataType, ShapeFunction, PSCorotatedLinear>;


/* Tetrahedral finite elements */
template<typename DataType>
using FEMPSNHTet = FEMPrincipalStretchTet<DataType, EnergyPSCoRot>; //Change EnergyPSCoRot to any other energy defined above to try out other marterials

typedef PhysicalSystemFEM<double, FEMPSNHTet> FEMLinearTets;

typedef World<double, std::tuple<FEMLinearTets *>,
std::tuple<ForceSpringFEMParticle<double> *, ForceParticlesGravity<double> *>,
std::tuple<ConstraintFixedPoint<double> *> > MyWorld;
typedef TimeStepperEulerImplicitLinear<double, AssemblerParallel<double, AssemblerEigenSparseMatrix<double> >,
AssemblerParallel<double, AssemblerEigenVector<double> > > MyTimeStepper;

typedef Scene<MyWorld, MyTimeStepper> MyScene;

void preStepCallback(MyWorld &world) {
    // This is an example callback
}


//Test newton solver by solving for static equilibrium of a bendy bar
int main(int argc, char **argv) {
    std::cout<<"Test Neohookean FEM \n";
    
    //Setup Physics
    MyWorld world;
    
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    
    readTetgen(V, F, dataDir()+"/meshesTetgen/Beam/Beam.node", dataDir()+"/meshesTetgen/Beam/Beam.ele");
    
    FEMLinearTets *test = new FEMLinearTets(V,F);
    world.addSystem(test);
    fixDisplacementMin(world, test);
    world.finalize(); //After this all we're ready to go (clean up the interface a bit later)
    
    auto q = mapStateEigen(world);
    q.setZero();
    
    Eigen::VectorXi indices = minVertices(test, 0);
    
    Eigen::SparseMatrix<double> P = fixedPointProjectionMatrix(indices, *test,world);
    
    MyTimeStepper stepper(0.05);
    
    //Display
    QGuiApplication app(argc, argv);
    
    MyScene *scene = new MyScene(&world, &stepper, preStepCallback);
    GAUSSVIEW(scene);
    
    std::cout << "Done." << std::endl;
    return app.exec();
    

}
