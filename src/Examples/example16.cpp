#include <functional>

#include <Qt3DIncludes.h>
#include <GaussIncludes.h>
#include <FEMIncludes.h>

//Any extra things I need such as constraints
#include <ConstraintFixedPoint.h>
#include <TimeStepperEulerImplicitLinear.h>
#include <Embedding.h>

using namespace Gauss;
using namespace FEM;
using namespace ParticleSystem; //For Force Spring

/* Tetrahedral finite elements */

//typedef physical entities I need

//typedef scene
typedef PhysicalSystemFEM<double, NeohookeanTet> FEMLinearTets;

typedef World<double, std::tuple< Embeddings::PhysicalSystemEmbeddedMesh<double, FEMLinearTets> *, PhysicalSystemParticleSingle<double> *>,
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
    
    //surface mesh
    Eigen::MatrixXd Vs;
    Eigen::MatrixXi Fs;
    
    
    readTetgen(V, F, dataDir()+"/meshesTetgen/Gargoyle/gargNested2.node", dataDir()+"/meshesTetgen/Gargoyle/gargNested2.ele");
    
    //Read mesh to embed
    if(!igl::readOBJ(dataDir()+"/meshes/OBJ/garg.obj", Vs, Fs))
        std::cout<<"Failed to read embedded mesh \n";
    
    //Embeddings::EmbeddingFunction<double, FEMLinearTets> embeddingFunction;
    Embeddings::PhysicalSystemEmbeddedMesh<double, FEMLinearTets> *test = new Embeddings::PhysicalSystemEmbeddedMesh<double, FEMLinearTets>(Vs,Fs, V, F);

    world.addSystem(test);
    fixDisplacementMin(world, test);
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
