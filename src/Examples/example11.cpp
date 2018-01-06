#include <functional>

#include <Qt3DIncludes.h>
#include <GaussIncludes.h>
#include <ParticleSystemIncludes.h>
#include <FEMIncludes.h>
#include <Newton.h>

#include <tuple>
//Any extra things I need such as constraints
#include <ConstraintFixedPoint.h>
#include <TimeStepperEulerImplicit.h>
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
typedef TimeStepperEulerImplicit<double, AssemblerEigenSparseMatrix<double>,
AssemblerEigenVector<double> > MyTimeStepper;

typedef Scene<MyWorld, MyTimeStepper> MyScene;


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
    
    readTetgen(V, F, dataDir()+"/meshesTetgen/Beam/Beam.node", dataDir()+"/meshesTetgen/Beam/Beam.ele");
    
    FEMLinearTets *test = new FEMLinearTets(V,F);
    world.addSystem(test);
    
    fixDisplacementMin(world, test);
    world.finalize(); //After this all we're ready to go (clean up the interface a bit later)
    
    auto q = mapStateEigen(world);
    q.setZero();
    
    MyTimeStepper stepper(0.1);
    
    //static solve minimizing the potential energy of an object
    
    //Newtons method requires energy, gradient, hessian, solver and initial position
    //Need some assemblers
    AssemblerEigenSparseMatrix<double> K;
    AssemblerEigenVector<double> f;
    SolverPardiso<Eigen::SparseMatrix<double, Eigen::RowMajor>> solver;
    
    //we're going to build equality constraints into our gradient and hessiang calcuations
    auto E = [&world](auto &a) { return getEnergy(world); }; //Not used by Newtons solver currently
    
    auto H = [&world, &K](auto &a)->auto & {
        //get mass matrix
        ASSEMBLEMATINIT(K, world.getNumQDotDOFs()+world.getNumConstraints(), world.getNumQDotDOFs()+world.getNumConstraints());
        ASSEMBLELIST(K, world.getSystemList(), getStiffnessMatrix);
        ASSEMBLELIST(K, world.getForceList(), getStiffnessMatrix);
        ASSEMBLELISTOFFSET(K, world.getConstraintList(), getGradient, world.getNumQDotDOFs(), 0);
        ASSEMBLELISTOFFSETTRANSPOSE(K, world.getConstraintList(), getGradient, 0, world.getNumQDotDOFs());
        ASSEMBLEEND(K);
    
        (*K) *= -1;
        return (*K);
    };
    
    auto g = [&world, &f](auto &a) -> auto & {
        ASSEMBLEVECINIT(f, world.getNumQDotDOFs()+world.getNumConstraints());
        ASSEMBLELIST(f, world.getForceList(), getForce);
        ASSEMBLELIST(f, world.getSystemList(), getForce);
        ASSEMBLEEND(f);
    
        (*f) *= -1;
        return (*f);
    };
    
     auto solve = [&solver](auto &A, auto &b)->auto & {
        return solver.solve(A,b);
    };
    
    auto update = [&world](auto &update) {
        mapStateEigen<0>(world) = update.head(world.getNumQDOFs());
    };
    
    //solve this using newton's method
    //auto x0 = Eigen::VectorXd(world.getNumQDotDOFs()+world.getNumConstraints());
    //x0.head(world.getNumQDotDOFs()) = mapStateEigen<0>(world);
    //Optimization::newton(E, g, H, solve, x0, update, 1e-5, 1000);
    
    //Display
    QGuiApplication app(argc, argv);
    
    MyScene *scene = new MyScene(&world, &stepper, preStepCallback);
    GAUSSVIEW(scene);
    
    return app.exec();
    
}
