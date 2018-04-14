#include <iostream>

#include <functional>

#include <Qt3DIncludes.h>
#include <GaussIncludes.h>
#include <FEMIncludes.h>

//Any extra things I need such as constraints
#include <ConstraintFixedPoint.h>
#include <CollisionDetector.h>
#include <CollisionsFCL.h>
#include <TimeStepperEulerImplicitLinearCollisions.h>
#include <Embedding.h>

//FCL
#include "fcl/math/bv/utility.h"
#include "fcl/narrowphase/collision.h"
#include "fcl/narrowphase/detail/gjk_solver_indep.h"
#include "fcl/narrowphase/detail/gjk_solver_libccd.h"
#include "fcl/narrowphase/detail/traversal/collision_node.h"

//Setup namespaces
using namespace Gauss;
using namespace FEM;
using namespace ParticleSystem; //For Force Spring
using namespace Collisions;

//Define convenient types
typedef ConstraintCollisionDetector<double, CollisionFCLImpl> MyCollisionDetector;
typedef PhysicalSystemFEM<double, NeohookeanTet> FEMNeohookeanTets;

typedef World<double, std::tuple< Embeddings::PhysicalSystemEmbeddedMesh<double, FEMNeohookeanTets> *>,
std::tuple<ForceSpringFEMParticle<double> *>,
std::tuple<ConstraintFixedPoint<double> *,MyCollisionDetector *> > MyWorld;

typedef TimeStepperEulerImplicitLinearCollisions<double, AssemblerEigenSparseMatrix<double>,
AssemblerEigenVector<double> > MyTimeStepper;

typedef Scene<MyWorld, MyTimeStepper> MyScene;

void preStepCallback(MyWorld &world) {
    // This is an example callback
}

int main(int argc, char **argv) {

    //load two tet meshes
    Eigen::MatrixXd V, Vs;
    Eigen::MatrixXi F, Fs;
    
    readTetgen(V, F, dataDir()+"/meshesTetgen/Beam/Beam.node", dataDir()+"/meshesTetgen/Beam/Beam.ele");
    
    //Extract boundary triangles
    Vs = V;
    
    //get the boundary facets for my data then write everything to disk
    igl::boundary_facets(F, Fs);
    Fs = Fs.rowwise().reverse().eval(); //igl boundary facets returns facet indices in reverse order
    
    //Embed
    Embeddings::PhysicalSystemEmbeddedMesh<double, FEMNeohookeanTets> *embeddedFEM = new Embeddings::PhysicalSystemEmbeddedMesh<double, FEMNeohookeanTets>(Vs,Fs, V, F);
    MyWorld world;

    world.addSystem(embeddedFEM);
    MyCollisionDetector cd(std::ref(world)); // this has to happen here, but I need to get rid of the circular dependence between collision detectors and the world
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
