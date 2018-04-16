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

typedef TimeStepperEulerImplicitLinearCollisions<double, AssemblerParallel<double, AssemblerEigenSparseMatrix<double> >,
AssemblerParallel<double, AssemblerEigenVector<double> > > MyTimeStepper;

typedef Scene<MyWorld, MyTimeStepper> MyScene;

void preStepCallback(MyWorld &world) {
    // This is an example callback
}

int main(int argc, char **argv) {

    //load two tet meshes
    Eigen::MatrixXd V, V2, Vs;
    Eigen::MatrixXi F, F2, Fs, Fs2;
    
    //readTetgen(V, F, dataDir()+"/meshesTetgen/Beam/Beam.node", dataDir()+"/meshesTetgen/Beam/Beam.ele");
    readTetgen(V, F, dataDir()+"/meshesTetgen/Gargoyle/gargNested2.node", dataDir()+"/meshesTetgen/Gargoyle/gargNested2.ele");
    readTetgen(V2, F2, dataDir()+"/meshesTetgen/Gargoyle/gargNested2.node", dataDir()+"/meshesTetgen/Gargoyle/gargNested2.ele");
    //readTetgen(V2, F2, dataDir()+"/meshesTetgen/Beam/Beam.node", dataDir()+"/meshesTetgen/Beam/Beam.ele");
    
    //Extract boundary triangles
    Vs = V;
    F2 = F;
    
    //get the boundary facets for my data then write everything to disk
    igl::boundary_facets(F, Fs);
    Fs = Fs.rowwise().reverse().eval(); //igl boundary facets returns facet indices in reverse order
    
    igl::boundary_facets(F2, Fs2);
    Fs2 = Fs2.rowwise().reverse().eval(); //igl boundary facets returns facet indices in reverse order
    
    V2.col(0) = (V2.col(0) + 1*Eigen::VectorXd::Ones(V2.rows())).eval(); //translate to the right by 5 m
    
    //Embed
    Embeddings::PhysicalSystemEmbeddedMesh<double, FEMNeohookeanTets> *embeddedFEM0 = new Embeddings::PhysicalSystemEmbeddedMesh<double, FEMNeohookeanTets>(V,Fs, V, F);
    Embeddings::PhysicalSystemEmbeddedMesh<double, FEMNeohookeanTets> *embeddedFEM1 = new Embeddings::PhysicalSystemEmbeddedMesh<double, FEMNeohookeanTets>(V2,Fs2, V2, F2);
    MyWorld world;

    world.addSystem(embeddedFEM0);
    world.addSystem(embeddedFEM1);
    MyCollisionDetector cd(std::ref(world)); // this has to happen here, but I need to get rid of the circular dependence between collision detectors and the world
    world.addInequalityConstraint(&cd);
    world.finalize();
    
    
    auto q = mapStateEigen(world);
    q.setZero();
    
    auto qDot = mapDOFEigen(embeddedFEM0->getQDot(), world);
    for(unsigned int ii=0; ii<qDot.rows(); ii+=3) qDot[ii] = 5.0;
    
    auto qDot2 = mapDOFEigen(embeddedFEM1->getQDot(), world);
    for(unsigned int ii=0; ii<qDot2.rows(); ii+=3) qDot2[ii] = -5.0;
    
    MyTimeStepper stepper(0.001);
    
    //Display
    QGuiApplication app(argc, argv);
    
    MyScene *scene = new MyScene(&world, &stepper, preStepCallback);
    GAUSSVIEW(scene);
    
    return app.exec();

}
