#include <Qt3DIncludes.h>
#include <GaussIncludes.h>
#include <FEMIncludes.h>

//Hexahedral FEM Test
#include <ConstraintFixedPoint.h>
#include <TimeStepperEulerImplicitLinear.h>

using namespace Gauss;
using namespace FEM;
using namespace ParticleSystem; //For Force Spring

/*Hexahedral finite elements */

//typedef physical entities I need

//typedef scene
typedef PhysicalSystemFEM<double, NeohookeanHex> FEMLinearHexes;

typedef World<double, std::tuple<FEMLinearHexes *>, std::tuple<ForceSpring<double> *>, std::tuple<ConstraintFixedPoint<double> *> > MyWorld;
typedef TimeStepperEulerImplictLinear<double, AssemblerEigenSparseMatrix<double>, AssemblerEigenVector<double> > MyTimeStepper;

typedef Scene<MyWorld, MyTimeStepper> MyScene;


int main(int argc, char **argv) {
    std::cout<<"Test neohookean material\n";
   
    Eigen::setNbThreads(1);
    
    std::cout<<"Test Linear FEM \n";
    
    //Setup Physics
    MyWorld world;
    
    
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    
    //Voxel grid from libigl
    igl::grid(Eigen::RowVector3i(50, 8, 8),  V);
    
    elementsFromGrid(Eigen::RowVector3i(50, 8,8), V, F);
    
    FEMLinearHexes *test = new FEMLinearHexes(V,F);
    
    world.addSystem(test);
    fixDisplacementMin(world, test);
    world.finalize(); //After this all we're ready to go (clean up the interface a bit later)
    
    auto q = mapStateEigen(world);
    q.setZero();
    
    AssemblerEigenSparseMatrix<double> massMatrix;
    AssemblerEigenVector<double> rhs;
    
    MyTimeStepper stepper(0.01);
    
    //Display
    QGuiApplication app(argc, argv);
    
    MyScene *scene = new MyScene(&world, &stepper);
    
    GAUSSVIEW(scene);
    gaussView.startScene();
    
    return app.exec();
    
    
    
}
