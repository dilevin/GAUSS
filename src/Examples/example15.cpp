#include <functional>

#include <Qt3DIncludes.h>
#include <GaussIncludes.h>
#include <ParticleSystemIncludes.h>
#include <FEMIncludes.h>

#include <tuple>

//Any extra things I need such as constraints
#include <ConstraintFixedPoint.h>
#include <TimeStepperEulerImplicitBFGS.h>
#include <type_traits>
using namespace Gauss;
using namespace FEM;
using namespace ParticleSystem; //For Force Spring

#include <type_traits>
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


void preStepCallback(MyWorld &world) {
    // This is an example callback
}


//Test newton solver by solving for static equilibrium of a bendy bar
int main(int argc, char **argv) {
    std::cout<<"Test Neohookean FEM \n";
    
    //Setup Physics
    MyWorld world;
    
    //new code -- load tetgen files
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    
    V.resize(4,3);
    F.resize(1,4);
    
    V << 0,0,0,
         1,0,0,
         0,1,0,
         0,0,1;
    F<<0,1,2,3;
    
    //simple test for my shape function matrix code
    FEMLinearTets *test = new FEMLinearTets(V,F);
    
    world.addSystem(test);
    
    world.finalize();
    
    Eigen::MatrixXd x;
    Eigen::VectorXi elements;
    x.resize(1,3);
    
    x<< 0.5, 0.5, 0.0;
    
    AssemblerEigenSparseMatrix<double> N;
    getShapeFunctionMatrix(N, elements, x, test->getImpl());
    
    std::cout<<"Shape function matrix: \n"<<(*N)<<"\n";
    return 0;
}
