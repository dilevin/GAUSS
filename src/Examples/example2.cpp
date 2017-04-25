#include <Qt3DIncludes.h>
#include <GaussIncludes.h>
#include <FEMIncludes.h>

//Any extra things I need such as constraints
#include <ConstraintFixedPoint.h>
#include <TimeStepperEulerImplicitLinear.h>

using namespace Gauss;
using namespace FEM;
using namespace ParticleSystem; //For Force Spring

//typedef physical entities I need

//typedef scene
typedef PhysicalSystemFEM<double, LinearTet> FEMLinearTets;
typedef World<double, std::tuple<FEMLinearTets *>, std::tuple<ForceSpring<double> *>, std::tuple<ConstraintFixedPoint<double> *> > MyWorld;
typedef TimeStepperEulerImplictLinear<double, AssemblerEigenSparseMatrix<double>,
AssemblerEigenVector<double> > MyTimeStepper;

typedef Scene<MyWorld, MyTimeStepper> MyScene;

//Utility functions to fix a bunch of points
template<typename World, typename FEMSystem>
void fixDisplacementLeft(World &world, FEMSystem *system) {
    //find all vertices with minimum x coordinate and fix DOF associated with them
    auto minX = system->getImpl().getV()(0,0);
    std::vector<unsigned int> minV;
    
    for(unsigned int ii=0; ii<system->getImpl().getV().rows(); ++ii) {
        
        if(system->getImpl().getV()(ii,0) < minX) {
            minX = system->getImpl().getV()(ii,0);
            minV.clear();
            minV.push_back(ii);
        } else if(fabs(system->getImpl().getV()(ii,0) - minX) < 1e-5) {
            minV.push_back(ii);
        }
    }
    
    //add a bunch of constraints
    for(auto iV : minV) {
        world.addConstraint(new ConstraintFixedPoint<decltype(minX)>(&system->getQ()[iV], Eigen::Vector3d(0,0,0)));
    }
}

int main(int argc, char **argv) {
    std::cout<<"Test Linear FEM \n";
    
    //Setup Physics
    MyWorld world;
    
    //new code -- load tetgen files
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    
    readTetgen(V, F, "../../data/meshesTetgen/Beam/Beam.node", "../../data/meshesTetgen/Beam/Beam.ele");
    
    std::cout<<"V:\n"<<V<<"\n\n";
    std::cout<<"F:\n"<<F<<"\n\n";
    //Mesh (old test code)
    //Eigen::MatrixXd V(4,3);
    //Eigen::MatrixXi F(1,4);
    
    //Load bar mesh and tetrahedralize using libigl
    //set boundary conditions and simulate using Gauss
    //V << 0, 0, 0,
    //1, 0, 0,
    //0, 1, 0,
    //0, 0, 1,
    
    //F<<0, 1, 2, 3;
    
    FEMLinearTets *test = new FEMLinearTets(V,F);
    
    //add my constraint
    //ConstraintFixedPoint<double> *fixedPoint = new ConstraintFixedPoint<double>(&test->getQ()[0], Eigen::Vector3d(0,0,0));
    //ConstraintFixedPoint<double> *fixedPoint1 = new ConstraintFixedPoint<double>(&test->getQ()[1], Eigen::Vector3d(0,0,0));
    //ConstraintFixedPoint<double> *fixedPoint2 = new ConstraintFixedPoint<double>(&test->getQ()[3], Eigen::Vector3d(0,0,0));
    world.addSystem(test);
    fixDisplacementLeft(world, test);
    //world.addConstraint(fixedPoint);
    //world.addConstraint(fixedPoint1);
    //world.addConstraint(fixedPoint2);
    world.finalize(); //After this all we're ready to go (clean up the interface a bit later)
    
    auto q = mapStateEigen(world);
    q.setZero();
    
    MyTimeStepper stepper(0.01);
    MyScene *scene = new MyScene(&world, &stepper);
    
    //Display
    QGuiApplication app(argc, argv);
    
    GAUSSVIEW(scene);
    gaussView.startScene();
    
    return app.exec();
}
