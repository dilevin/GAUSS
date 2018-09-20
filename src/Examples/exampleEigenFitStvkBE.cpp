#include <functional>

#include <Qt3DIncludes.h>
#include <GaussIncludes.h>
#include <FEMIncludes.h>

//Any extra things I need such as constraints
#include <ConstraintFixedPoint.h>
#include <TimeStepperEigenFitStvkSMWBE.h>
#include <EigenFitStvk.h>
#include <fstream>
#include <igl/boundary_facets.h>
#include <igl/volume.h>
#include <igl/dihedral_angles.h>

using namespace Gauss;
using namespace FEM;
using namespace ParticleSystem; //For Force Spring

/* Tetrahedral finite elements */

//typedef physical entities I need

//typedef scene
//typedef PhysicalSystemFEM<double, StvkHFixedTet> FEMLinearTets;
typedef PhysicalSystemFEM<double, StvkHFixedTet> FEMLinearTets;

typedef World<double, std::tuple<FEMLinearTets *>,
std::tuple<ForceSpringFEMParticle<double> *, ForceParticlesGravity<double> *>,
std::tuple<ConstraintFixedPoint<double> *> > MyWorld;

//typedef World<double, std::tuple<FEMLinearTets *,PhysicalSystemParticleSingle<double> *>,
//                      std::tuple<ForceSpringFEMParticle<double> *>,
//                      std::tuple<ConstraintFixedPoint<double> *> > MyWorld;
//typedef TimeStepperEigenFitStvkSMW<double, AssemblerEigenSparseMatrix<double>, AssemblerEigenVector<double>> MyTimeStepper;
typedef TimeStepperEigenFitStvkSMWBE<double, AssemblerParallel<double, AssemblerEigenSparseMatrix<double>>, AssemblerParallel<double, AssemblerEigenVector<double>> > MyTimeStepper;

typedef Scene<MyWorld, MyTimeStepper> MyScene;


//typedef TimeStepperEigenFitStvkSI<double, AssemblerParallel<double, AssemblerEigenSparseMatrix<double> >,
//AssemblerParallel<double, AssemblerEigenVector<double> > > MyTimeStepper;

//typedef Scene<MyWorld, MyTimeStepper> MyScene;

// used for preStepCallback. should be delete
std::vector<ConstraintFixedPoint<double> *> movingConstraints;
std::vector<ConstraintFixedPoint<double> *> fixedConstraints;
Eigen::VectorXi movingVerts;
Eigen::VectorXi fixedVerts;
Eigen::MatrixXd V, Vtemp;
Eigen::MatrixXi F;
Eigen::MatrixXi surfF;
Eigen::MatrixXi surfFf;

char **arg_list;
unsigned int istep;

void preStepCallback(MyWorld &world) {
}

int main(int argc, char **argv) {
    //    arg list
    //    1: full path to coarse mesh
    //    2: full path to fine mesh
    //    3: youngs modulus (SI unit)
    //    4: constraint threshold (for defualt constraint profile)
    //    5: constraint profile switch
    //    6: name of the initial deformation profile
    //    7. number of time steps
    //    8. flag for using hausdorff distance
    //    9. number of modes to modifies
    //    10. constraint direction
    //    11. step size
    
    
    std::cout<<"Test Stvk FEM EigenFitStvk\n";
    
    //Setup Physics
    MyWorld world;
    
    arg_list = argv;
    
    Eigen::MatrixXd Vf;
    Eigen::MatrixXi Ff;
    
    
    //    define the file separator base on system
    const char kPathSeparator =
#ifdef _WIN32
    '\\';
#else
    '/';
#endif
    
    if (argc > 1) {
        // must supply all 9 parameters
        
        std::string cmeshname = argv[1];
        std::string fmeshname = argv[2];
        
        readTetgen(V, F, dataDir()+cmeshname+".node", dataDir()+cmeshname+".ele");
        readTetgen(Vf, Ff, dataDir()+fmeshname+".node", dataDir()+fmeshname+".ele");
        
        Vtemp = V;
        //        find the surface mesh
        igl::boundary_facets(F,surfF);
        igl::boundary_facets(Ff,surfFf);
        
        std::string::size_type found = cmeshname.find_last_of(kPathSeparator);
        //    acutal name for the mesh, no path
        std::string cmeshnameActual = cmeshname.substr(found+1);
        
        //    acutal name for the mesh, no path
        std::string fmeshnameActual = fmeshname.substr(found+1);
        
        
        //    parameters
        double youngs = atof(argv[3]);
        double poisson = 0.45;
        double constraint_tol = atof(argv[4]);
        int const_profile = atoi(argv[5]);
        int initial_def = atoi(argv[6]);
        int numSteps = atoi(argv[7]);
        bool hausdorff = atoi(argv[8]);
        int numModes = atoi(argv[9]);
        int constraint_dir = atoi(argv[10]); // constraint direction. 0 for x, 1 for y, 2 for z
        double step_size = atof(argv[11]);
        bool dynamic_flag = atoi(argv[12]);
        //
        // send the constraint switch in as well, or the fine embedded mesh. ugly
        // the flag indicate whether to recalculated or not
        // need to pass the material and constraint parameters to embedding too. need to do it again below. ugly
        // also use the last two args to determine how many modes to fix. have to put it here now. ugly
        EigenFitStvk *test = new EigenFitStvk(V,F,Vf,Ff,dynamic_flag,youngs,poisson,constraint_dir,constraint_tol, const_profile,hausdorff,numModes,cmeshnameActual,fmeshnameActual);
        
        
        world.addSystem(test);
        
        
        // projection matrix for constraints
        Eigen::SparseMatrix<double> P;
        
        // constraint profile should be 2
        if(const_profile == 2)
        {
            //            zero gravity
            Eigen::Vector3x<double> g;
            g(0) = 0;
            g(1) = 0;
            g(2) = 0;
            
            for(unsigned int iel=0; iel<test->getImpl().getF().rows(); ++iel) {
                
                test->getImpl().getElement(iel)->setGravity(g);
                
            }
            
            Eigen::loadMarketVector(movingVerts, "def_init/" + cmeshnameActual + "_fixed_min_verts.mtx");
            
            for(unsigned int ii=0; ii<movingVerts.rows(); ++ii) {
                movingConstraints.push_back(new ConstraintFixedPoint<double>(&test->getQ()[movingVerts[ii]], Eigen::Vector3d(0,0,0)));
                world.addConstraint(movingConstraints[ii]);
            }
            fixDisplacementMin(world, test,constraint_dir,constraint_tol);
            
            world.finalize(); //After this all we're ready to go (clean up the interface a bit later)
            
            P = fixedPointProjectionMatrix(movingVerts, *test,world);
        }
        else if (const_profile == 4 || const_profile == 5 || const_profile == 6 || const_profile == 7 || const_profile == 8)
        {
            //            zero gravity
            Eigen::Vector3x<double> g;
            g(0) = 0;
            g(1) = 0;
            g(2) = 0;
            
            for(unsigned int iel=0; iel<test->getImpl().getF().rows(); ++iel) {
                
                test->getImpl().getElement(iel)->setGravity(g);
                
            }
            // read constraints
            
            //            Eigen::loadMarketVector(fixedVerts, "def_init/" + cmeshnameActual + "_fixed_min_verts.mtx");
            //
            movingVerts = minVertices(test, constraint_dir, constraint_tol);//indices for moving parts
            
            for(unsigned int ii=0; ii<movingVerts.rows(); ++ii) {
                movingConstraints.push_back(new ConstraintFixedPoint<double>(&test->getQ()[movingVerts[ii]], Eigen::Vector3d(0,0,0)));
                world.addConstraint(movingConstraints[ii]);
            }
            fixDisplacementMin(world, test,constraint_dir,constraint_tol);
            
            world.finalize(); //After this all we're ready to go (clean up the interface a bit later)
            
            P = fixedPointProjectionMatrix(movingVerts, *test,world);
        }
        else
        {
            std::cout<<"warning: wrong constraint profile\n";
        }
        
        
        // set material
        for(unsigned int iel=0; iel<test->getImpl().getF().rows(); ++iel) {
            
            test->getImpl().getElement(iel)->setParameters(youngs, poisson);
            
        }
        
        // initialize the state (position and velocity)
        auto q = mapStateEigen(world);
        //        std::cout<<q.rows()<<std::endl;
        //        std::cout<<q.cols()<<std::endl;
        //        q.setZero();
        
        
        if (strcmp(argv[6],"0")==0) {
            // if specified no initial deformation
            q.setZero();
        }
        else
        {
            std::cout<<"warning: wrong initial deformation\n";
        }
        
        MyTimeStepper stepper(step_size,P,numModes);
        
        //         the number of steps to take
        
        unsigned int file_ind = 0;
        unsigned int mode = 0;
        unsigned int vert_idx = 0;
        //        Eigen::MatrixXd coarse_eig_def;
        //        Eigen::VectorXd fine_eig_def;
        
        struct stat buf;
        unsigned int idxc;
        
        for(istep=0; istep<numSteps ; ++istep) {
            try{
                stepper.step(world);
                
                throw test->flag;
            }
            catch(int e)
            {
                if(e == 2)
                {
                    std::cout<<"hausdorff distance check failed\n";
                    return 1;
                }
            }
            // acts like the "callback" block for moving constraint
            if (const_profile == 2)
            {
                // constraint profile 2 will move some vertices
                //script some motion
                for(unsigned int jj=0; jj<movingConstraints.size(); ++jj) {
                    
                    auto v_q = mapDOFEigen(movingConstraints[jj]->getDOF(0), world.getState());
                    Eigen::Vector3d new_q = (istep)*Eigen::Vector3d(0.0,-1.0/100,0.0);
                    v_q = new_q;
                    
                }
            }
            else if (const_profile == 4 )
            {
                for(unsigned int jj=0; jj<movingConstraints.size(); ++jj) {
                    
                    auto v_q = mapDOFEigen(movingConstraints[jj]->getDOF(0), world.getState());
                    //
                    if ((istep) < 50) {
                        Eigen::Vector3d new_q = (istep)*Eigen::Vector3d(-1.0/100,0.0,0.0);
                        v_q = new_q;
                    }
                    
                }
            }
            else if (const_profile == 5)
            {
                for(unsigned int jj=0; jj<movingConstraints.size(); ++jj) {
                    
                    auto v_q = mapDOFEigen(movingConstraints[jj]->getDOF(0), world.getState());
                    //
                    if ((istep) < 50) {
                        Eigen::Vector3d new_q = (istep)*Eigen::Vector3d(0.0,-1.0/100,0.0);
                        v_q = new_q;
                    }
                    
                }
            }else if (const_profile == 6)
            {
                for(unsigned int jj=0; jj<movingConstraints.size(); ++jj) {
                    
                    auto v_q = mapDOFEigen(movingConstraints[jj]->getDOF(0), world.getState());
                    //
                    if ((istep) < 50) {
                        Eigen::Vector3d new_q = (istep)*Eigen::Vector3d(0.0,0.0,-1.0/100);
                        v_q = new_q;
                    }
                    
                    
                }
            }else if (const_profile == 7)
            {
                Eigen::VectorXd Xvel;
                Eigen::loadMarketVector(Xvel, "mouseXvel.mtx");
                Eigen::VectorXd Yvel;
                Eigen::loadMarketVector(Yvel, "mouseYvel.mtx");
                Eigen::VectorXd Zvel;
                Eigen::loadMarketVector(Zvel, "mouseZvel.mtx");
                
                for(unsigned int jj=0; jj<movingConstraints.size(); ++jj) {
                    
                    auto v_q = mapDOFEigen(movingConstraints[jj]->getDOF(0), world.getState());
                    //
                    if ((istep) < 250) {
                        //                        Eigen::Vector3d new_q = (istep)*Eigen::Vector3d(0.0,0.0,-1.0/100);
                        v_q(0) += 0.1*Xvel(istep);
                        v_q(1) += 0.1*Yvel(istep);
                        v_q(2) += 0.1*Zvel(istep);
                    }
                    
                    
                }
            }else if (const_profile == 8)
            {
                Eigen::VectorXd Xvel;
                Eigen::loadMarketVector(Xvel, "mouseXvel.mtx");
                Eigen::VectorXd Yvel;
                Eigen::loadMarketVector(Yvel, "mouseYvel.mtx");
                Eigen::VectorXd Zvel;
                Eigen::loadMarketVector(Zvel, "mouseZvel.mtx");
                
                for(unsigned int jj=0; jj<movingConstraints.size(); ++jj) {
                    
                    auto v_q = mapDOFEigen(movingConstraints[jj]->getDOF(0), world.getState());
                    //
                    if ((istep) < 250) {
                        //                        Eigen::Vector3d new_q = (istep)*Eigen::Vector3d(0.0,0.0,-1.0/100);
                        if(Xvel(istep) <= 0){   v_q(0) += std::max(Xvel(istep),-0.005);}
                        else{ v_q(0) += std::min(Xvel(istep),0.005);}
                        if(Yvel(istep) <= 0){   v_q(1) += std::max(Yvel(istep),-0.005);}
                        else{ v_q(0) += std::min(Yvel(istep),0.005);}
                        if(Zvel(istep) <= 0){   v_q(2) += std::max(Zvel(istep),-0.005);}
                        else{ v_q(0) += std::min(Zvel(istep),0.005);}
                        //
                        //                            v_q(1) += Yvel(istep);
                        //                        v_q(2) += Zvel(istep);
                    }
                    
                    
                }
            }
            //output data stream into text
            // the following ones append one number to an opened file along the simulation
            std::ofstream ofile;
            ofile.open("KE.txt", std::ios::app); //app is append which means it will put the text at the end
            ofile << std::get<0>(world.getSystemList().getStorage())[0]->getImpl().getKineticEnergy(world.getState()) << std::endl;
            ofile.close();
            
            //output data stream into text
            ofile.open("PE.txt", std::ios::app); //app is append which means it will put the text at the end
            ofile << std::get<0>(world.getSystemList().getStorage())[0]->getImpl().getStrainEnergy(world.getState()) << std::endl;
            ofile.close();
            
            //output data stream into text
            ofile.open("Hamiltonian.txt", std::ios::app); //app is append which means it will put the text at the end
            ofile << std::get<0>(world.getSystemList().getStorage())[0]->getImpl().getEnergy(world.getState()) << std::endl;
            ofile.close();
            
            // check if the file already exist
            std::string filename = "pos" + std::to_string(file_ind) + ".obj";
            while (stat(filename.c_str(), &buf) != -1)
            {
                file_ind++;
                filename = "pos" + std::to_string(file_ind) + ".obj";
            }
            
            // rest pos for the coarse mesh getGeometry().first is V
            q = mapStateEigen(world);
            idxc = 0;
            Eigen::MatrixXd V_disp = std::get<0>(world.getSystemList().getStorage())[0]->getGeometry().first;
            // output mesh position with only surface mesh
            //            igl::writeOBJ("surfpos_rest" + std::to_string(file_ind) + ".obj",V_disp,surfF);
            
            
            // get the mesh position
            for(unsigned int vertexId=0;  vertexId < std::get<0>(world.getSystemList().getStorage())[0]->getGeometry().first.rows(); ++vertexId) {
                
                V_disp(vertexId,0) += q(idxc);
                idxc++;
                V_disp(vertexId,1) += q(idxc);
                idxc++;
                V_disp(vertexId,2) += q(idxc);
                idxc++;
            }
            
            // output mesh position with elements
            igl::writeOBJ("pos" + std::to_string(file_ind) + ".obj",V_disp,std::get<0>(world.getSystemList().getStorage())[0]->getGeometry().second);
            
            // output mesh position with only surface mesh
            igl::writeOBJ("surfpos" + std::to_string(file_ind) + ".obj",V_disp,surfF);
            
            if (numModes != 0) {
                // declare variable for fine mesh rest pos
                Eigen::MatrixXd Vf_disp;
                
                // embedded V
                auto fine_q = mapStateEigen<0>(test->getFineWorld());
                //                std::cout<<"N rows: " <<(*(test->N)).rows() <<"\n";
                //                std::cout<<"N cols: " <<(*(test->N)).cols() <<"\n";
                //                std::cout<<"q rows: " <<q.rows() <<"\n";
                //                std::cout<<"fine q rows: " <<fine_q.rows() <<"\n";
                fine_q = (*(test->N)) * q.head(q.rows()/2);
                idxc = 0; // reset index counter
                Vf_disp = std::get<0>(test->getFineWorld().getSystemList().getStorage())[0]->getGeometry().first;
                // output mesh position with only surface mesh
                for(unsigned int vertexId=0;  vertexId < std::get<0>(test->getFineWorld().getSystemList().getStorage())[0]->getGeometry().first.rows(); ++vertexId) {
                    
                    Vf_disp(vertexId,0) += fine_q(idxc);
                    idxc++;
                    Vf_disp(vertexId,1) += fine_q(idxc);
                    idxc++;
                    Vf_disp(vertexId,2) += fine_q(idxc);
                    idxc++;
                }
                // output mesh position with only surface mesh
                igl::writeOBJ("finesurfpos" + std::to_string(file_ind) + ".obj",Vf_disp,surfFf);
                
            }
        }
        
    }
    else
    {
        // using all default paramters for EigenFitStvk
        
        //    default example meshes
        std::string cmeshname = "/meshesTetWild/brick/brick_surf_4";
        std::string fmeshname = "/meshesTetWild/brick/brick_surf_4";
        
        readTetgen(V, F, dataDir()+cmeshname+".node", dataDir()+cmeshname+".ele");
        readTetgen(Vf, Ff, dataDir()+fmeshname+".node", dataDir()+fmeshname+".ele");
        
        std::string::size_type found = cmeshname.find_last_of(kPathSeparator);
        //    acutal name for the mesh, no path
        std::string cmeshnameActual = cmeshname.substr(found+1);
        
        
        //    default parameters
        double youngs = 2e5;
        double poisson = 0.45;
        int constraint_dir = 0; // constraint direction. 0 for x, 1 for y, 2 for z
        double constraint_tol = 1e-2;
        bool dynamic_flag = true;
        int const_profile = 0;
        bool hausdorff = false;
        int numModes = 3;
        
        // no constraint switch so just create the EigenFitStvk obj with constraint switch set to 1
        // the flag indicate whether to recalculated or not
        // need to pass the material and constraint parameters to embedding too. need to do it again below. ugly
        // also use the last two args to determine how many modes to fix. default not using hausdorff distance, and use 10 modes. have to put it here now.  ugly
        //        EigenFitStvk *test = new EigenFitStvk(V,F,Vf,Ff,dynamic_flag,youngs,poisson,constraint_dir,constraint_tol, const_profile,hausdorff,numModes,cmeshnameActual,fmeshnameActual);
        
        EigenFitStvk *test = new EigenFitStvk(V,F,Vf,Ff,true,youngs,poisson,constraint_dir,constraint_tol, const_profile,hausdorff,numModes," "," ");
        
        // set material
        for(unsigned int iel=0; iel<test->getImpl().getF().rows(); ++iel) {
            
            test->getImpl().getElement(iel)->setParameters(youngs, poisson);
            
        }
        
        world.addSystem(test);
        
        world.finalize(); //After this all we're ready to go (clean up the interface a bit later)
        // IMPORTANT, need to finalized before fix boundary
        
        //    default constraint
        fixDisplacementMin(world, test,constraint_dir,constraint_tol);
        
        // construct the projection matrix for stepper
        Eigen::VectorXi indices = minVertices(test, constraint_dir,constraint_tol);
        Eigen::SparseMatrix<double> P = fixedPointProjectionMatrix(indices, *test,world);
        
        
        auto q = mapStateEigen(world);
        
        //    default to zero deformation
        q.setZero();
        
        
        MyTimeStepper stepper(0.01,P,0);
        
        //Display
        QGuiApplication app(argc, argv);
        
        MyScene *scene = new MyScene(&world, &stepper, preStepCallback);
        GAUSSVIEW(scene);
        
        return app.exec();
        
    }
    
}
