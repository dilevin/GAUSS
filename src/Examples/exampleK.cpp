#include <Qt3DIncludes.h>
#include <GaussIncludes.h>
#include <FEMIncludes.h>
#include <igl/writeOBJ.h>
#include <igl/readOBJ.h>

//Plame strain triangular element test
#include <ConstraintFixedPoint.h>
#include <TimeStepperEulerImplicitLinear.h>

using namespace Gauss;
using namespace FEM;
using namespace ParticleSystem; //For Force Spring

//typedef scene
typedef PhysicalSystemFEM<double, LinearPlaneStrainTri> FEMPlaneStrainTri;

typedef World<double, std::tuple<FEMPlaneStrainTri *>, std::tuple<ForceSpringFEMParticle<double> *>, std::tuple<ConstraintFixedPoint<double> *> > MyWorld;
typedef TimeStepperEulerImplicitLinear<double, AssemblerEigenSparseMatrix<double>, AssemblerEigenVector<double> > MyTimeStepper;

typedef Scene<MyWorld, MyTimeStepper> MyScene;
MyTimeStepper stepper(1e-3);
void preStepCallback(MyWorld &world) {
    // stepper.step(world);
    // std::cout<<"hi"<<std::endl;
}
int main(int argc, char **argv) {
    std::cout<<"Test neohookean material\n";

    Eigen::setNbThreads(1);

    std::cout<<"Test Linear FEM \n";

    //Setup Physics
    MyWorld world;


    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    // V.resize(9,3);
    // V<< 0, 0, 0,
    // 8, 0, 0,
    // 8, 8, 0,
    // 4, 0, 0,
    // 8, 4, 0,
    // 4, 4, 0,
    // 0, 8, 0,
    // 4, 8, 0,
    // 0, 4, 0;
    // F.resize(8, 3);
    // F<<0, 4, 8,
    // 0, 8, 7,
    // 4, 1, 5,
    // 4, 8, 8,
    // 8, 5, 2,
    // 8, 2, 6,
    // 7, 8, 6,
    // 7, 6, 3;
    //simple square subdivided into two triangles
    igl::readOBJ("/home/vismay/Basis/mesh1.obj", V, F);
    std::cout<<"V F"<<std::endl;
    std::cout<<V<<std::endl;
    std::cout<<F<<std::endl;
    FEMPlaneStrainTri *test = new FEMPlaneStrainTri(V,F);
    test->getImpl().getElement(0)->setParameters(1e6, 0.45);


    world.addSystem(test);
    fixDisplacementMin(world, test);
    world.finalize(); //After this all we're ready to go (clean up the interface a bit later)

    auto q = mapStateEigen(world);
    q.setZero();
    auto p = mapDOFEigen(test->getQ(), world);
    AssemblerParallel<double, AssemblerEigenVector<double>> f;
    getForceVector(f, world);
    std::cout<<*f<<std::endl;
    std::cout<<p<<std::endl;
    AssemblerParallel<double, AssemblerEigenSparseMatrix<double> > stiffness;
    AssemblerParallel<double, AssemblerEigenSparseMatrix<double> > mass;
    getMassMatrix(mass, world);
    getStiffnessMatrix(stiffness, world);
    std::cout<<"Mass"<<std::endl;
    std::cout<<*mass<<std::endl;
    std::cout<<"Stiffness"<<std::endl;
    std::cout<<*stiffness<<std::endl;


    // for(int i=0; i<1000; ++i){
    //   stepper.step(world);
    //   if(i%100==0)
    //   {
    //     auto p = mapDOFEigen(test->getQ(), world);
    //     Eigen::Map<Eigen::MatrixXd> dV(p.data(), V.rows(), V.cols());
    //     Eigen::MatrixXd nV = V+dV;
    //     igl::writeOBJ("mesh2_"+std::to_string(i)+".obj", nV, F);
    //   }
    // }

    //Display
    QGuiApplication app(argc, argv);

    MyScene *scene = new MyScene(&world, &stepper, preStepCallback);

    GAUSSVIEW(scene);
    gaussView.startScene();
    return app.exec();



}
