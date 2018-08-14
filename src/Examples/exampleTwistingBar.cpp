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
#include <igl/boundary_facets.h>
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
typedef TimeStepperEulerImplicit<double, AssemblerParallel<double, AssemblerEigenSparseMatrix<double>>,
AssemblerParallel<double, AssemblerEigenVector<double> > > MyTimeStepper;

typedef Scene<MyWorld, MyTimeStepper> MyScene;


std::vector<ConstraintFixedPoint<double> *> movingConstraints;
Eigen::MatrixXd V, Vtemp;
Eigen::MatrixXi F;
Eigen::VectorXi movingVerts;
Eigen::MatrixXi surfF;

int current_frame = 0;

void preStepCallback(MyWorld &world) {
    current_frame++;
    // This is an example callback

    //script some motion
    Eigen::AngleAxis<double> rot(0.01 * current_frame, Eigen::Vector3d(1.0,0.0,0.0));

    for(unsigned int jj=0; jj<movingConstraints.size(); ++jj) {
        //if(movingConstraints[jj]->getImpl().getFixedPoint()[0] > -3) {
            Eigen::Vector3d v = V.row(movingVerts[jj]);
            // rot * v is the new position, so new_u is new velocity direction
            Eigen::Vector3d new_u = rot * v - v - mapDOFEigen(movingConstraints[jj]->getDOF(0), world.getState());
            movingConstraints[jj]->getImpl().setFixedPoint(rot*v, new_u/0.1);
        //}
    }
    auto q = mapStateEigen(world);
    Eigen::saveMarketVector(q, "brick_surf_2_def" + std::to_string(current_frame) + ".mtx");
    
    
    unsigned int idxc = 0;
    
    // get the mesh position
    for(unsigned int vertexId=0;  vertexId < std::get<0>(world.getSystemList().getStorage())[0]->getGeometry().first.rows(); ++vertexId) {
        
        Vtemp(vertexId,0) += q(idxc);
        idxc++;
        Vtemp(vertexId,1) += q(idxc);
        idxc++;
        Vtemp(vertexId,2) += q(idxc);
        idxc++;
    }
    
    // output mesh position with elements
    igl::writeOBJ("savedpos" + std::to_string(current_frame) +".obj",Vtemp,surfF);
    igl::writeOBJ("savedpos0" + std::to_string(current_frame) +".obj",V,surfF);
//    Eigen::saveMarketVector(q,"saveddef.mtx");

}

int main(int argc, char **argv) {
    std::cout<<"Test Neohookean FEM \n";

    //Setup Physics
    MyWorld world;

    //new code -- load tetgen files
    readTetgen(V, F, dataDir()+"/meshesTetWild/brick/brick_surf_2.node", dataDir()+"/meshesTetWild/brick/brick_surf_2.ele");
    Vtemp = V;
    
    igl::boundary_facets(F,surfF);
    
    FEMLinearTets *test = new FEMLinearTets(V,F);

    world.addSystem(test);
    
    fixDisplacementMin(world, test,0,1e-2); //fix one side
    Eigen::saveMarketVector(minVertices(test,0,1e-2),"brick_surf_2_fixed_min_verts.mtx");
    movingVerts = maxVertices(test, 0,1e-1);//indices for moving parts

    for(unsigned int ii=0; ii<movingVerts.rows(); ++ii) {
        // initialize the moving constraint by setting the position. velocity is initialized to zero, but will be changed in callback
        movingConstraints.push_back(new ConstraintFixedPoint<double>(&test->getQ()[movingVerts[ii]], Eigen::Vector3d(0,0,0)));
        world.addConstraint(movingConstraints[ii]);
    }

    world.finalize(); //After this all we're ready to go (clean up the interface a bit later)

    auto q = mapStateEigen(world);
    q.setZero();

    MyTimeStepper stepper(0.1, 200);

    //Display
    QGuiApplication app(argc, argv);

    MyScene *scene = new MyScene(&world, &stepper, preStepCallback);
    GAUSSVIEW(scene);

    return app.exec();

}
