#include <functional>

#include <Qt3DIncludes.h>
#include <GaussIncludes.h>
#include <FEMIncludes.h>

//Any extra things I need such as constraints
#include <ConstraintFixedPoint.h>
#include <TimeStepperEigenFitSMW.h>
#include <EigenFit.h>

using namespace Gauss;
using namespace FEM;
using namespace ParticleSystem; //For Force Spring

/* Tetrahedral finite elements */

//typedef physical entities I need

//typedef scene
typedef PhysicalSystemFEM<double, NeohookeanHFixedTet> FEMLinearTets;

typedef World<double, std::tuple<FEMLinearTets *>,
std::tuple<ForceSpringFEMParticle<double> *, ForceParticlesGravity<double> *>,
std::tuple<ConstraintFixedPoint<double> *> > MyWorld;

//typedef World<double, std::tuple<FEMLinearTets *,PhysicalSystemParticleSingle<double> *>,
//                      std::tuple<ForceSpringFEMParticle<double> *>,
//                      std::tuple<ConstraintFixedPoint<double> *> > MyWorld;
//typedef TimeStepperEigenFitSMW<double, AssemblerEigenSparseMatrix<double>, AssemblerEigenVector<double>> MyTimeStepper;
typedef TimeStepperEigenFitSMW<double, AssemblerParallel<double, AssemblerEigenSparseMatrix<double>>, AssemblerParallel<double, AssemblerEigenVector<double>> > MyTimeStepper;

typedef Scene<MyWorld, MyTimeStepper> MyScene;


//typedef TimeStepperEigenFitSI<double, AssemblerParallel<double, AssemblerEigenSparseMatrix<double> >,
//AssemblerParallel<double, AssemblerEigenVector<double> > > MyTimeStepper;

//typedef Scene<MyWorld, MyTimeStepper> MyScene;


void preStepCallback(MyWorld &world) {
    // This is an example callback
}

int main(int argc, char **argv) {
    std::cout<<"Test Neohookean FEM EigenFit\n";
    
    //Setup Physics
    MyWorld world;
    
    //new code -- load tetgen files
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    
    Eigen::MatrixXd Vf;
    Eigen::MatrixXi Ff;
    
    std::string cmeshname = "/meshesTetgen/arma/arma_6";
    std::string fmeshname = "/meshesTetgen/arma/arma_1";
    
    if (argc > 1) {
        cmeshname = argv[1];
        fmeshname = argv[2];
    }
    
    readTetgen(V, F, dataDir()+cmeshname+".node", dataDir()+cmeshname+".ele");

    readTetgen(Vf, Ff, dataDir()+fmeshname+".node", dataDir()+fmeshname+".ele");

//    readTetgen(Vf, Ff, dataDir()+"/meshesTetgen/bar/barmesh8.node", dataDir()+"/meshesTetgen/bar/barmesh8.ele");
//    readTetgen(Vf, Ff, dataDir()+"/meshesTetgen/Beam/Beam.node", dataDir()+"/meshesTetgen/Beam/Beam.ele");

    double youngs = 5e5;
    double poisson = 0.45;
    int constraint_dir = 2; // constraint direction. 0 for x, 1 for y, 2 for z
    double constraint_tol = 2e-1;
    
    // the flag indicate whether to recalculated or not
    // need to pass the material and constraint parameters to embedding too. need to do it again below. ugly
    EigenFit *test = new EigenFit(V,F,Vf,Ff,true,youngs,poisson,constraint_dir,constraint_tol);
    
    for(unsigned int iel=0; iel<test->getImpl().getF().rows(); ++iel) {
        
        test->getImpl().getElement(iel)->setParameters(youngs, poisson);
        
    }
//    test->getImpl.getElements.setParameters(1e5,0.45);
    world.addSystem(test);
    fixDisplacementMin(world, test,constraint_dir,constraint_tol);
    world.finalize(); //After this all we're ready to go (clean up the interface a bit later)
    
     auto q = mapStateEigen(world);
     q.setZero();
    
    
    Eigen::VectorXi indices = minVertices(test, constraint_dir,constraint_tol);
    Eigen::SparseMatrix<double> P = fixedPointProjectionMatrix(indices, *test,world);
    
    
     MyTimeStepper stepper(0.01,P);
    if(argc > 3) {
        
        unsigned int file_ind = 0;
        std::string name = "pos";
        std::string fformat = ".obj";
        std::string filename = name + std::to_string(file_ind) + fformat;
        struct stat buf;
        unsigned int idx;
        
        for(unsigned int istep=0; istep<atoi(argv[3]) ; ++istep) {
            stepper.step(world);
            
            //output data here
            std::ofstream ofile;
            ofile.open("KE.txt", std::ios::app); //app is append which means it will put the text at the end
            ofile << std::get<0>(world.getSystemList().getStorage())[0]->getImpl().getKineticEnergy(world.getState()) << std::endl;
            ofile.close();
            
            while (stat(filename.c_str(), &buf) != -1)
            {
                file_ind++;
                filename = name + std::to_string(file_ind) + fformat;
            }
            
                idx = 0;
                // getGeometry().first is V
                Eigen::MatrixXd V_disp = std::get<0>(world.getSystemList().getStorage())[0]->getGeometry().first;
            
                for(unsigned int vertexId=0;  vertexId < std::get<0>(world.getSystemList().getStorage())[0]->getGeometry().first.rows(); ++vertexId) {
            
                    // because getFinePosition is in EigenFit, not another physical system Impl, so don't need getImpl()
                    V_disp(vertexId,0) += q(idx);
                    idx++;
                    V_disp(vertexId,1) += q(idx);
                    idx++;
                    V_disp(vertexId,2) += q(idx);
                    idx++;
                }
                igl::writeOBJ(filename,V_disp,std::get<0>(world.getSystemList().getStorage())[0]->getGeometry().second);

        }
    }
    else {
    //Display
    QGuiApplication app(argc, argv);
    
     MyScene *scene = new MyScene(&world, &stepper, preStepCallback);
     GAUSSVIEW(scene);
    
    return app.exec();
    }
}
