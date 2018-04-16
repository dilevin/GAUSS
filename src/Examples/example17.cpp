#include <functional>
#include <Qt3DIncludes.h>
#include <GaussIncludes.h>
#include <ParticleSystemIncludes.h>
#include <FEMIncludes.h>
#include <CollisionDetector.h>
#include <CollisionsFloor.h>

#include <tuple>
//Any extra things I need such as constraints
#include <ConstraintFixedPoint.h>
#include <TimeStepperEulerImplicitLinearCollisions.h>
#include <type_traits>
using namespace Gauss;
using namespace Collisions;
using namespace FEM;
using namespace ParticleSystem;

typedef PhysicalSystemFEM<double, NeohookeanTet> FEMLinearTets;
typedef ConstraintCollisionDetector<double, CollisionFloorImpl> MyCollisionDetector;

//World using embedded mesh FEM
typedef World<double, std::tuple< Embeddings::PhysicalSystemEmbeddedMesh<double, FEMLinearTets> *, PhysicalSystemParticleSingle<double> *>,
std::tuple<ForceSpringFEMParticle<double> *, ForcePoint<double> *>,
std::tuple<ConstraintFixedPoint<double> *, MyCollisionDetector *> > MyWorld;

typedef TimeStepperEulerImplicitLinearCollisions<double, AssemblerEigenSparseMatrix<double>,
AssemblerEigenVector<double> > MyTimeStepper;

typedef Scene<MyWorld, MyTimeStepper> MyScene;

//Per vertex force
Eigen::VectorXd perVertexForce;

void preStepCallback(MyWorld &world) {
    // This is an example callback
}

int main(int argc, char **argv) {
    std::cout<<"Collisions Test \n";
    
    MyWorld world;
    
    //Storage for vertices and facets for simulation mesh
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
    
    MyCollisionDetector cd(std::ref(world), Eigen::Vector3d(-1,1.0,0.0), Eigen::Vector3d(0.0, -5,0.0));
    world.addSystem(test);
    world.addInequalityConstraint(&cd);
    world.finalize();
    
    auto q = mapStateEigen(world);
    q.setZero();
    
    MyTimeStepper stepper(0.01);
    
    //Display
    QGuiApplication app(argc, argv);
    
    MyScene *scene = new MyScene(&world, &stepper, preStepCallback);
    GAUSSVIEW(scene);
    
    return app.exec();
    
    
}
