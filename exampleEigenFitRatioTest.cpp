#include <functional>

#include <Qt3DIncludes.h>
#include <GaussIncludes.h>
#include <FEMIncludes.h>

//Any extra things I need such as constraints
#include <ConstraintFixedPoint.h>
#include <StepperEigenFitRatioTest.h>
#include <EigenFit.h>
#include <fstream>
#include <igl/boundary_facets.h>

using namespace Gauss;
using namespace FEM;
using namespace ParticleSystem; //For Force Spring

/* Tetrahedral finite elements */

//typedef physical entities I need

//typedef scene
//typedef PhysicalSystemFEM<double, NeohookeanHFixedTet> FEMLinearTets;
typedef PhysicalSystemFEM<double, NeohookeanHFixedTet> FEMLinearTets;

typedef World<double, std::tuple<FEMLinearTets *>,
std::tuple<ForceSpringFEMParticle<double> *, ForceParticlesGravity<double> *>,
std::tuple<ConstraintFixedPoint<double> *> > MyWorld;

//typedef World<double, std::tuple<FEMLinearTets *,PhysicalSystemParticleSingle<double> *>,
//                      std::tuple<ForceSpringFEMParticle<double> *>,
//                      std::tuple<ConstraintFixedPoint<double> *> > MyWorld;
//typedef TimeStepperEigenFitSMW<double, AssemblerEigenSparseMatrix<double>, AssemblerEigenVector<double>> MyTimeStepper;
typedef StepperEigenFitRatioTest<double, AssemblerParallel<double, AssemblerEigenSparseMatrix<double>>, AssemblerParallel<double, AssemblerEigenVector<double>> > MyTimeStepper;

typedef Scene<MyWorld, MyTimeStepper> MyScene;


//typedef TimeStepperEigenFitSI<double, AssemblerParallel<double, AssemblerEigenSparseMatrix<double> >,
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
    //    // This is an example callback
    //    if (atoi(arg_list[5]) == 2)
    //    {
    //        // This is an example callback
    //
    //        //script some motion
    ////
    //        if (istep < 50) {
    //
    //            for(unsigned int jj=0; jj<movingConstraints.size(); ++jj) {
    //                if(movingConstraints[jj]->getImpl().getFixedPoint()[0] > -3) {
    //                    Eigen::Vector3d v = V.row(movingVerts[jj]);
    //                    Eigen::Vector3d new_p = v + Eigen::Vector3d(0.0,1.0/10,0.0);
    //                    movingConstraints[jj]->getImpl().setFixedPoint(new_p);
    //                }
    //            }
    //        }
    //    }
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
    
    
    std::cout<<"Test Neohookean FEM EigenFit\n";
    
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
    int constraint_dir = atoi(argv[10]); // constraint direction. 0 for x, 1 for y, 2 for z
    double constraint_tol = atof(argv[4]);
    int numModes = atoi(argv[9]);
    bool hausdorff = atoi(argv[8]);
    int const_profile = atoi(argv[5]);
    int initial_def = atoi(argv[6]);
    double step_size = atof(argv[11]);
    int numSteps = atoi(argv[7]);
    bool dynamic_flag = atoi(argv[12]);
    //
    // send the constraint switch in as well, or the fine embedded mesh. ugly
    // the flag indicate whether to recalculated or not
    // need to pass the material and constraint parameters to embedding too. need to do it again below. ugly
    // also use the last two args to determine how many modes to fix. have to put it here now. ugly
    EigenFit *test = new EigenFit(V,F,Vf,Ff,dynamic_flag,youngs,poisson,constraint_dir,constraint_tol, const_profile,hausdorff,numModes,cmeshnameActual,fmeshnameActual);
    
    
    world.addSystem(test);
    
    
    // projection matrix for constraints
    Eigen::SparseMatrix<double> P;
    
    if ((const_profile) == 0) {
        // constraint switch
        
        //            zero gravity
        Eigen::Vector3x<double> g;
        g(0) = 0;
        g(1) = 0;
        g(2) = 0;
        
        for(unsigned int iel=0; iel<test->getImpl().getF().rows(); ++iel) {
            
            test->getImpl().getElement(iel)->setGravity(g);
            
        }
        
        world.finalize(); //After this all we're ready to go (clean up the interface a bit later)
        
        //            set the projection matrix to identity because there is no constraint to project
        //            Eigen::SparseMatrix<double> P;
        P.resize(V.rows()*3,V.rows()*3);
        P.setIdentity();
        //            std::cout<<P.rows();
        //            no constraints
    }
    else if(const_profile == 1)
    {
        //    default constraint
        fixDisplacementMin(world, test,constraint_dir,constraint_tol);
        world.finalize(); //After this all we're ready to go (clean up the interface a bit later)
        
        // construct the projection matrix for stepper
        Eigen::VectorXi indices = minVertices(test, constraint_dir,constraint_tol);
        P = fixedPointProjectionMatrix(indices, *test,world);
        // for debugging
        //            if(numModes!=0)
        //            {
        //                Eigen::saveMarket(P, "coarse_P.dat");
        //            }
        //            else
        //            {
        //                Eigen::saveMarket(P, "noeigenfit_P.dat");
        //            }
    }
    else if (const_profile == 2)
    {
        
        
        movingVerts = minVertices(test, constraint_dir, constraint_tol);//indices for moving parts
        //
        for(unsigned int ii=0; ii<movingVerts.rows(); ++ii) {
            movingConstraints.push_back(new ConstraintFixedPoint<double>(&test->getQ()[movingVerts[ii]], Eigen::Vector3d(0,0,0)));
            world.addConstraint(movingConstraints[ii]);
        }
        fixDisplacementMin(world, test,constraint_dir,constraint_tol);
        
        world.finalize(); //After this all we're ready to go (clean up the interface a bit later)
        
        P = fixedPointProjectionMatrix(movingVerts, *test,world);
        
    }
    else if (const_profile == 3)
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
        
        Eigen::loadMarketVector(fixedVerts, "def_init/" + cmeshnameActual + "_fixed_min_verts.mtx");
        //
        for(unsigned int ii=0; ii<fixedVerts.rows(); ++ii) {
            fixedConstraints.push_back(new ConstraintFixedPoint<double>(&test->getQ()[fixedVerts[ii]], Eigen::Vector3d(0,0,0)));
            world.addConstraint(fixedConstraints[ii]);
        }
        
        world.finalize(); //After this all we're ready to go (clean up the interface a bit later)
        std::cout<<"# of fixedVerts" << fixedVerts.rows()<<std::endl;
        P = fixedPointProjectionMatrix(fixedVerts, *test,world);
        
    }
    
    // set material
    for(unsigned int iel=0; iel<test->getImpl().getF().rows(); ++iel) {
        
        test->getImpl().getElement(iel)->setParameters(youngs, poisson);
        
    }
    
    // initialize the state (position and velocity)
    auto q = mapStateEigen(world);
    //        if (numModes != 0) {
    //        Eigen::VectorXd fine_q(Vf.rows() * 3);
    //        std::cout<<mapStateEigen(test->getFineWorld())<<std::endl;
    //        fine_q = mapStateEigen(test->getFineWorld());
    //        }
    //    default to zero deformation
    std::cout<<q.rows()<<std::endl;
    std::cout<<q.cols()<<std::endl;
    //        q.setZero();
    
    
    
    MyTimeStepper stepper(step_size,P,numModes,cmeshnameActual,fmeshnameActual);
    
    //         the number of steps to take
    
    unsigned int file_ind = 0;
    unsigned int mode = 0;
    unsigned int vert_idx = 0;
    //        Eigen::MatrixXd coarse_eig_def;
    //        Eigen::VectorXd fine_eig_def;
    
    struct stat buf;
    unsigned int idxc;
    
    for(istep=0; istep<numSteps ; ++istep) {
        stepper.step(world);
        
        // acts like the "callback" block for moving constraint
        if (const_profile == 2)
        {
            // constraint profile 2 will move some vertices
            //script some motion
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
        
    }
}




