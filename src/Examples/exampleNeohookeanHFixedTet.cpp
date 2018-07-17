#include <functional>

#include <Qt3DIncludes.h>
#include <GaussIncludes.h>
#include <FEMIncludes.h>

//Any extra things I need such as constraints
#include <ConstraintFixedPoint.h>
#include <TimeStepperEulerImplicitLinearCProjection.h>
#include <sys/stat.h>

using namespace Gauss;
using namespace FEM;
using namespace ParticleSystem; //For Force Spring

/* Tetrahedral finite elements */

//typedef physical entities I need

//typedef scene
typedef PhysicalSystemFEM<double, NeohookeanHFixedTet> FEMLinearTets;

//typedef World<double, std::tuple<FEMLinearTets *>,
//std::tuple<ForceSpringFEMParticle<double> *, ForceParticlesGravity<double> *>,
//std::tuple<ConstraintFixedPoint<double> *> > MyWorld;

//typedef PhysicalSystemFEM<double, NeohookeanHFixedTet> FEMLinearTets;


typedef World<double, std::tuple<FEMLinearTets *,PhysicalSystemParticleSingle<double> *>,
std::tuple<ForceSpringFEMParticle<double> *>,
std::tuple<ConstraintFixedPoint<double> *> > MyWorld;
typedef TimeStepperEulerImplicitLinearCProjection<double, AssemblerParallel<double, AssemblerEigenSparseMatrix<double>>, AssemblerParallel<double, AssemblerEigenVector<double>>> MyTimeStepper;

typedef Scene<MyWorld, MyTimeStepper> MyScene;


// used for preStepCallback. should be delete
std::vector<ConstraintFixedPoint<double> *> movingConstraints;
Eigen::VectorXi movingVerts;
Eigen::MatrixXd V;
Eigen::MatrixXi F;
char **arg_list;
unsigned int istep;


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
    
    //    define the file separator base on system
    const char kPathSeparator =
#ifdef _WIN32
    '\\';
#else
    '/';
#endif
    
    if (argc > 1) {
        
        std::string meshname = argv[1];
        
        readTetgen(V, F, dataDir()+meshname+".node", dataDir()+meshname+".ele");
        
        std::string::size_type found = meshname.find_last_of(kPathSeparator);
        //    acutal name for the mesh, no path
        std::string meshnameActual = meshname.substr(found+1);
        
        //    parameters
        double youngs = atof(argv[2]);
        double poisson = 0.45;
        int constraint_dir = atoi(argv[8]); // constraint direction. 0 for x, 1 for y, 2 for z
        double constraint_tol = atof(argv[3]);
        
        FEMLinearTets *test = new FEMLinearTets(V,F);
        
        world.addSystem(test);
        
        // projection matrix for constraints
        Eigen::SparseMatrix<double> P;
        if (atoi(argv[4]) == 0) {
            //             constraint switch
            
            //            zero gravity
            Eigen::Vector3x<double> g;
            g(0) = 0;
            g(1) = 0;
            g(2) = 0;
            
            for(unsigned int iel=0; iel<test->getImpl().getF().rows(); ++iel) {
                
                test->getImpl().getElement(iel)->setGravity(g);
                
            }
            
            world.finalize(); //After this all we're ready to go (clean up the interface a bit later)
            P.resize(V.rows()*3,V.rows()*3);
            P.setIdentity();
            
        }
        else if(atoi(argv[4]) == 1)
        {
            //    default constraint
            fixDisplacementMin(world, test,constraint_dir,constraint_tol);
            //            world.finalize(); //After this all we're ready to go (clean up the interface a bit later)
            world.finalize(); //After this all we're ready to go (clean up the interface a bit later)
            
            // construct the projection matrix for stepper
            Eigen::VectorXi indices = minVertices(test, constraint_dir,constraint_tol);
            P = fixedPointProjectionMatrix(indices, *test,world);
            
            
        }
        else if (atoi(argv[4]) == 2)
        {
            //            //            zero gravity
            //            Eigen::Vector3x<double> g;
            //            g(0) = 0;
            //            g(1) = 0;
            //            g(2) = 0;
            //
            //            for(unsigned int iel=0; iel<test->getImpl().getF().rows(); ++iel) {
            //
            //                test->getImpl().getElement(iel)->setGravity(g);
            //
            //            }
            
            movingVerts = minVertices(test, constraint_dir, constraint_tol);//indices for moving parts
            
            for(unsigned int ii=0; ii<movingVerts.rows(); ++ii) {
                movingConstraints.push_back(new ConstraintFixedPoint<double>(&test->getQ()[movingVerts[ii]], Eigen::Vector3d(0,0,0)));
                world.addConstraint(movingConstraints[ii]);
            }
            fixDisplacementMin(world, test,constraint_dir,constraint_tol);
            
            world.finalize(); //After this all we're ready to go (clean up the interface a bit later)
            
            P = fixedPointProjectionMatrix(movingVerts, *test,world);
            
        }
        
        // set material
        for(unsigned int iel=0; iel<test->getImpl().getF().rows(); ++iel) {
            
            test->getImpl().getElement(iel)->setParameters(youngs, poisson);
            
        }
        
        auto q = mapStateEigen(world);
        
        //    default to zero deformation
        q.setZero();
        
        if (strcmp(argv[5],"0")==0) {
            
            q.setZero();
        }
        else
        {
            
            std::string qfileName(argv[5]);
            Eigen::VectorXd  tempv;
            loadMarketVector(tempv,qfileName);
            //            std::cout<<tempv.size();
            //            std::cout<<tempv;
            q = tempv;
        }
        
        
        MyTimeStepper stepper(0.01,P);
        
        //         the number of steps to take
        
        unsigned int file_ind = 0;
        std::string name = "pos";
        std::string fformat = ".obj";
        std::string filename = name + std::to_string(file_ind) + fformat;
        std::string qname = meshnameActual + "q";
        std::string qfformat = ".mtx";
        std::string qfilename = qname + std::to_string(file_ind) + qfformat;
        
        struct stat buf;
        unsigned int idx;
        
        for(unsigned int istep=0; istep<atoi(argv[6]) ; ++istep) {
            stepper.step(world);
            
            
            // acts like the "callback" block
            if (atoi(argv[4]) == 2)
            {
                // This is an example callback
                
                //script some motion
                //
                
                
                for(unsigned int jj=0; jj<movingConstraints.size(); ++jj) {
                    
                    auto v_q = mapDOFEigen(movingConstraints[jj]->getDOF(0), world.getState());
//
//                    if ((istep%150) < 50) {
//                        Eigen::Vector3d new_q = (istep%150)*Eigen::Vector3d(0.0,-1.0/100,0.0);
//                        v_q = new_q;
//                    }
//                    else if ((istep%150) < 100)
//                    {}
//                    else
//                    {
//                        Eigen::Vector3d new_q =  (150-(istep%150))*Eigen::Vector3d(0.0,-1.0/100,0.0);
//                        v_q = new_q;
//                    }
                    Eigen::Vector3d new_q = (istep)*Eigen::Vector3d(0.0,-1.0/100,0.0);
                    v_q = new_q;
                }
            }
            
            
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
            //            saveMarketVector(q, qfilename);
        }
        world.addSystem(test);
        
        
    }
    else
    {
        // using all default paramters for eigenfit
        
        //    default example meshes
        std::string meshname = "/meshesTetgen/arma/arma_6";
        
        readTetgen(V, F, dataDir()+meshname+".node", dataDir()+meshname+".ele");
        
        std::string::size_type found = meshname.find_last_of(kPathSeparator);
        //    acutal name for the mesh, no path
        std::string meshnameActual = meshname.substr(found+1);
        
        
        //    default parameters
        double youngs = 5e5;
        double poisson = 0.45;
        int constraint_dir = atoi(argv[8]); // constraint direction. 0 for x, 1 for y, 2 for z
        double constraint_tol = 2e-1;
        
        FEMLinearTets *test = new FEMLinearTets(V,F);
        
        // set material
        for(unsigned int iel=0; iel<test->getImpl().getF().rows(); ++iel) {
            
            test->getImpl().getElement(iel)->setParameters(youngs, poisson);
            
        }
        
        world.addSystem(test);
        //    default constraint
        fixDisplacementMin(world, test,constraint_dir,constraint_tol);
        world.finalize(); //After this all we're ready to go (clean up the interface a bit later)
        // IMPORTANT, need to finalized before fix boundary
        //
        //        //    default constraint
        //        fixDisplacementMin(world, test,constraint_dir,constraint_tol);
        
        // construct the projection matrix for stepper
        Eigen::VectorXi indices = minVertices(test, constraint_dir,constraint_tol);
        Eigen::SparseMatrix<double> P = fixedPointProjectionMatrix(indices, *test,world);
        
        
        auto q = mapStateEigen(world);
        
        //    default to zero deformation
        q.setZero();
        
        
        MyTimeStepper stepper(0.01,P);
        
        //Display
        QGuiApplication app(argc, argv);
        
        MyScene *scene = new MyScene(&world, &stepper, preStepCallback);
        GAUSSVIEW(scene);
        
        return app.exec();
    }
    //
    //    std::string meshname = "/meshesTetgen/arma/arma_6";
    //
    //    if (argc > 1) {
    //        meshname = argv[1];
    //    }
    //
    //    readTetgen(V, F, dataDir()+meshname+".node", dataDir()+meshname+".ele");
    //
    //    FEMLinearTets *test = new FEMLinearTets(V,F);
    //
    //    double youngs = 5e5;
    //    double poisson = 0.45;
    //    int constraint_dir = 0; // constraint direction. 0 for x, 1 for y, 2 for z
    //    double constraint_tol = 2e-1;
    //
    //    if (argc > 2) {
    //        youngs = atof(argv[2]);
    //    }
    //
    //    if (argc > 3) {
    //        constraint_tol = atof(argv[3]);
    //    }
    //
    //    for(unsigned int iel=0; iel<test->getImpl().getF().rows(); ++iel) {
    //
    //        test->getImpl().getElement(iel)->setParameters(youngs, poisson);
    //
    //    }
    //
    //    world.addSystem(test);
    //    fixDisplacementMin(world, test, constraint_dir,constraint_tol);
    //    world.finalize(); //After this all we're ready to go (clean up the interface a bit later)
    //
    //    auto q = mapStateEigen(world);
    //    q.setZero();
    //
    //    MyTimeStepper stepper(0.01);
    //    if(argc > 4) {
    //
    //        unsigned int file_ind = 0;
    //        std::string name = "pos";
    //        std::string fformat = ".obj";
    //        std::string filename = name + std::to_string(file_ind) + fformat;
    //        struct stat buf;
    //        unsigned int idx;
    //
    //        for(unsigned int istep=0; istep<atoi(argv[4]) ; ++istep) {
    //            stepper.step(world);
    //
    //            //output data here
    //            std::ofstream ofile;
    //            ofile.open("KE.txt", std::ios::app); //app is append which means it will put the text at the end
    //            ofile << std::get<0>(world.getSystemList().getStorage())[0]->getImpl().getKineticEnergy(world.getState()) << std::endl;
    //            ofile.close();
    //
    //            while (stat(filename.c_str(), &buf) != -1)
    //            {
    //                file_ind++;
    //                filename = name + std::to_string(file_ind) + fformat;
    //            }
    //
    //            idx = 0;
    //            // getGeometry().first is V
    //            Eigen::MatrixXd V_disp = std::get<0>(world.getSystemList().getStorage())[0]->getGeometry().first;
    //
    //            for(unsigned int vertexId=0;  vertexId < std::get<0>(world.getSystemList().getStorage())[0]->getGeometry().first.rows(); ++vertexId) {
    //
    //                // because getFinePosition is in EigenFit, not another physical system Impl, so don't need getImpl()
    //                V_disp(vertexId,0) += q(idx);
    //                idx++;
    //                V_disp(vertexId,1) += q(idx);
    //                idx++;
    //                V_disp(vertexId,2) += q(idx);
    //                idx++;
    //            }
    //            igl::writeOBJ(filename,V_disp,std::get<0>(world.getSystemList().getStorage())[0]->getGeometry().second);
    //
    //        }
    //    }
    //    else {
    //        //Display
    //        QGuiApplication app(argc, argv);
    //
    //        MyScene *scene = new MyScene(&world, &stepper, preStepCallback);
    //        GAUSSVIEW(scene);
    //
    //        return app.exec();
    //    }
}
