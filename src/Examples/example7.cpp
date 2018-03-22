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
typedef TimeStepperEulerImplicitLinear<double, AssemblerEigenSparseMatrix<double>, AssemblerEigenVector<double>  > MyTimeStepper;

typedef Scene<MyWorld, MyTimeStepper> MyScene;

Eigen::MatrixXi triangulateRect(Eigen::MatrixXd V) {
    int res = sqrt(V.rows());
    Eigen::MatrixXi F((res - 1)*(res-1) * 2, 3);
    
    int idx = 0;
    for (int r = 0; r < res-1; r++) {
        for (int c = 0; c < res-1; c++) {
            int rowFront = res*(r);
            int nextRowFront = res*(r + 1);
            
            {
                int t3 = nextRowFront + c;
                int t2 = rowFront + c + 1;
                int t1 = rowFront + c;
                F.row(idx) << t1, t2, t3;
                idx++;
            }
            {
                int t3 = nextRowFront + c;
                int t2 = nextRowFront + c + 1;
                int t1 = rowFront + c + 1;
                F.row(idx) << t1, t2, t3;
                idx++;
            }
            
        }
        
    }
    return F;
}

int main(int argc, char **argv) {
    std::cout<<"Test neohookean material\n";

    Eigen::setNbThreads(1);

    std::cout<<"Test Linear FEM \n";

    //Setup Physics
    MyWorld world;


    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::MatrixXd fieldNodes;
    
    igl::grid(Eigen::RowVector3i(100, 100, 1),  fieldNodes);
    
    V.resize(fieldNodes.rows(), fieldNodes.cols());
    V << fieldNodes;
    
    int res = sqrt(V.rows());
    F.resize((res - 1)*(res - 1) * 2, 3);
    F << triangulateRect(fieldNodes);
    
    // Add a third, Z column of zeros.
    V.conservativeResize(fieldNodes.rows(), 3);
    V.col(2).setZero();
    
    //simple square subdivided into two triangles
    /*V.resize(4,3);
    V << 0,0,0,
         1,0,0,
         1,1,0,
         0,1,0;

    F.resize(2,3);
    F << 0, 1, 2,
         0, 2, 3;*/

    FEMPlaneStrainTri *test = new FEMPlaneStrainTri(V,F);

    world.addSystem(test);
    fixDisplacementMin(world, test);
    world.finalize(); //After this all we're ready to go (clean up the interface a bit later)

    auto q = mapStateEigen(world);
    q.setZero();

    MyTimeStepper stepper(0.01, false);

    //Display
    QGuiApplication app(argc, argv);
    
    MyScene *scene = new MyScene(&world, &stepper);
    
    GAUSSVIEW(scene);
    
    return app.exec();


    return 1;



}
