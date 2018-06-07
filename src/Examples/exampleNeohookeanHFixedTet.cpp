#include <functional>

#include <Qt3DIncludes.h>
#include <GaussIncludes.h>
#include <FEMIncludes.h>

//Any extra things I need such as constraints
#include <ConstraintFixedPoint.h>
#include <TimeStepperEulerImplicitLinear.h>
#include <sys/stat.h>

using namespace Gauss;
using namespace FEM;
using namespace ParticleSystem; //For Force Spring

/* Tetrahedral finite elements */

//typedef physical entities I need

//typedef scene
typedef PhysicalSystemFEM<double, NeohookeanHFixedTet> FEMLinearTets;

typedef World<double, std::tuple<FEMLinearTets *,PhysicalSystemParticleSingle<double> *>,
                      std::tuple<ForceSpringFEMParticle<double> *>,
                      std::tuple<ConstraintFixedPoint<double> *> > MyWorld;
typedef TimeStepperEulerImplicitLinear<double, AssemblerParallel<double, AssemblerEigenSparseMatrix<double>>, AssemblerParallel<double, AssemblerEigenVector<double>>> MyTimeStepper;

typedef Scene<MyWorld, MyTimeStepper> MyScene;


void preStepCallback(MyWorld &world) {
    // This is an example callback
}

int main(int argc, char **argv) {
    std::cout<<"Test NeohookeanHFixed FEM \n";
    
    //Setup Physics
    MyWorld world;
    
    //new code -- load tetgen files
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    
    std::string meshname = "/meshesTetgen/arma/arma_6";
    
    if (argc > 1) {
        meshname = argv[1];
    }
    
    readTetgen(V, F, dataDir()+meshname+".node", dataDir()+meshname+".ele");
   
    FEMLinearTets *test = new FEMLinearTets(V,F);

    
    for(unsigned int iel=0; iel<test->getImpl().getF().rows(); ++iel) {
        
        test->getImpl().getElement(iel)->setParameters(5e5, 0.45);
        
    }

    world.addSystem(test);
    fixDisplacementMin(world, test, 2,2e-1);
    world.finalize(); //After this all we're ready to go (clean up the interface a bit later)
    
    auto q = mapStateEigen(world);
    q.setZero();
    
    MyTimeStepper stepper(0.01);
    if(argc > 2) {
        
        unsigned int file_ind = 0;
        std::string name = "pos";
        std::string fformat = ".obj";
        std::string filename = name + std::to_string(file_ind) + fformat;
        struct stat buf;
        unsigned int idx;
        
        for(unsigned int istep=0; istep<atoi(argv[2]) ; ++istep) {
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
