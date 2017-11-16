#include <Qt3DIncludes.h>
#include <GaussIncludes.h>
#include <FEMIncludes.h>

//Plame strain triangular element test
#include <ConstraintFixedPoint.h>
#include <TimeStepperEulerImplicitLinear.h>

using namespace Gauss;
using namespace FEM;
using namespace ParticleSystem; //For Force Spring

//typedef scene
typedef PhysicalSystemFEM<double, LinearPlaneStrainTri> FEMPlaneStrainTri;

typedef World<double, std::tuple<FEMPlaneStrainTri *>, std::tuple<ForceSpringFEMParticle<double> *>, std::tuple<ConstraintFixedPoint<double> *> > MyWorld;
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
    
    //simple square subdivided into two triangles
    V.resize(4,3);
    V << 0,0,0,
         1,0,0,
         1,1,0,
         0,1,0;
    
    F.resize(2,3);
    F << 0, 1, 2,
         0, 2, 3;
    
    FEMPlaneStrainTri *test = new FEMPlaneStrainTri(V,F);
    
    world.addSystem(test);
    fixDisplacementMin(world, test);
    world.finalize(); //After this all we're ready to go (clean up the interface a bit later)
    
    auto q = mapStateEigen(world);
    q.setZero();
    
    MyTimeStepper stepper(0.01);
    
    stepper.step(world);
    
    std::cout<<q<<"\n";
    
    
    return 1;
    
    
    
}
