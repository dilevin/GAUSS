#include "mex.h"
#include "class_handle.hpp"

//Gauss Stuff
#include <GaussIncludes.h>
#include <FEMIncludes.h>
#include <ForceSpring.h>
#include <ForceParticleGravity.h>

//Any extra things I need such as constraints
#include <ConstraintFixedPoint.h>
#include <TimeStepperEulerImplicitLinear.h>

//Utilities
#include "UtilitiesEigenMex.h"

using namespace Gauss;
using namespace FEM;
using namespace ParticleSystem; //For Force Spring

/* Tetrahedral finite elements */

//typedef physical entities I need

//typedef scene
typedef PhysicalSystemFEM<double, LinearTet> FEMLinearTets;
typedef PhysicalSystemFEM<double, LinearPlaneStrainTri> FEMPlaneStrainTri;
typedef PhysicalSystemFEM<double, NeohookeanTet> FEMNeohookeanTets;

//Generic world object that supports the sim objects I need
typedef World<  double,
                std::tuple<
                    FEMLinearTets *,
                    FEMPlaneStrainTri *,
                    FEMNeohookeanTets *
                >,
                std::tuple<ForceSpringFEMParticle<double> *, ForceParticlesGravity<double> *>,
                std::tuple<ConstraintFixedPoint<double> *>
        > WorldFEM;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    
    Eigen::setNbThreads(1); //Eigens open mp breaks things
    
    // Get the command string
    char cmd[64];
	if (nrhs < 1 || mxGetString(prhs[0], cmd, sizeof(cmd)))
		mexErrMsgTxt("First input should be a command string less than 64 characters long.");
 
    // New
    if (!strcmp("new", cmd)) {
        if(nrhs < 4)
            mexErrMsgTxt("Arguments: need FEM Type (string), Vertex (matrix) and Face array (matrix)");
    
        // Check parameters
        if (nlhs != 1)
            mexErrMsgTxt("New: One output expected.");
        
        // Return a handle to a new C++ instance
        
        WorldFEM *world = new WorldFEM();
        
        char * femType = mxArrayToString(prhs[1]);
        
        if(!femType) {
            mexPrintf("Invalid Parameter: FEM Type not a string. Please clear and reinitialize FEM object. \n");
            return;
        }
        
        mexPrintf("%s\n", femType);
        if(strcmp(femType, "elastic_linear_tetrahedra") == 0) {
            mexPrintf("Initialize Linear Elastic Tetrahedra\n");
            FEMLinearTets *FEM = new FEMLinearTets(matlabToDouble(prhs[2]), matlabToInt32(prhs[3]));
            world->addSystem(FEM);
        } else if(strcmp(femType, "elastic_linear_plane_strain_tri")==0) {
            mexPrintf("Initialize Linear Elastic Plane Strain Tri\n");
            FEMPlaneStrainTri *FEM = new FEMPlaneStrainTri(matlabToDouble(prhs[2]), matlabToInt32(prhs[3]));
            world->addSystem(FEM);
        } else if(strcmp(femType, "neohookean_linear_tetrahedra")==0) {
            mexPrintf("Initialize Neohookean Elastic Tetrahedra\n");
            FEMNeohookeanTets *FEM = new FEMNeohookeanTets(matlabToDouble(prhs[2]), matlabToInt32(prhs[3]));
            world->addSystem(FEM);
        } else {
           mexPrintf("Invalid Physical System. GAUSS not initalized \n");
        }
        
        world->finalize(); //After this all we're ready to go (clean up the interface a bit later)
        auto q = mapStateEigen(*world);
        q.setZero();
        plhs[0] = convertPtr2Mat<WorldFEM>(world);
        return;
    }
    
    // Check there is a second input, which should be the class instance handle
    if (nrhs < 2)
		mexErrMsgTxt("Second input should be a class instance handle.");
    
    // Delete
    if (strcmp("delete", cmd) == 0) {
        // Destroy the C++ object
        destroyObject<WorldFEM>(prhs[1]);
        // Warn if other commands were ignored
        if (nlhs != 0 || nrhs != 2)
            mexWarnMsgTxt("Delete: Unexpected arguments ignored.");
        return;
    }
    
    // Get the class instance pointer from the second input
    WorldFEM *dummy_instance = convertMat2Ptr<WorldFEM>(prhs[1]);
    
    // Call the various class methods
    // Train    
    if (!strcmp("state", cmd)) {
        // Check parameters
        if (nlhs < 0 || nrhs < 1)
            mexErrMsgTxt("Not enough parameters to state.");
        Eigen::Map<Eigen::VectorXd> state = mapStateEigen(*dummy_instance);
        
        //copy state into matlab vector
        mwSize dims[2];
        dims[0] = state.rows();
        dims[1] = 1;
        mxArray *returnData = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
        memcpy(mxGetPr(returnData), state.data(), sizeof(double)*state.rows());
        plhs[0] = returnData;
        return;
    }
    
    if (!strcmp("setState", cmd)) {
        // Check parameters
        if (nlhs < 0 || nrhs < 3)
            mexErrMsgTxt("Not enough parameters to state.");
        
        double* A = 0;
        size_t  m = 0;
        size_t  n = 0;
        
        A = mxGetPr(prhs[2]);
        // get the dimensions of the first parameter
        m = mxGetM(prhs[2]);
        n = mxGetN(prhs[2]);
        
        Eigen::Map<Eigen::VectorXd> state = mapStateEigen(*dummy_instance);
        
        if(m != state.rows()) {
            mexPrintf("Input vector is of size %i x %i but state is of size %i \n", m,n, state.rows(), 1);
        }
        
        state = Eigen::Map<Eigen::VectorXd>(A, state.rows());
        
        //copy state into matlab vector
        //mwSize dims[2];
        //dims[0] = state.rows();
        //dims[1] = 1;
        //mxArray *returnData = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
        //memcpy(mxGetPr(returnData), state.data(), sizeof(double)*state.rows());
        //plhs[0] = returnData;
        return;
    }
    
    //set position level DOFs
    if (!strcmp("setQ", cmd)) {
        
        // Check parameters
        if (nlhs < 0 || nrhs < 3)
            mexErrMsgTxt("Not enough parameters to state.");
        
        double* A = 0;
        size_t  m = 0;
        size_t  n = 0;
        
        A = mxGetPr(prhs[2]);
        // get the dimensions of the first parameter
        m = mxGetM(prhs[2]);
        n = mxGetN(prhs[2]);
        
        Eigen::Map<Eigen::VectorXd> state = mapStateEigen<0>(*dummy_instance);
        
        if(m != state.rows()) {
            mexPrintf("Input vector is of size %i x %i but state is of size %i \n", m,n, state.rows(), 1);
        }
        
        state = Eigen::Map<Eigen::VectorXd>(A, state.rows());
        
        //copy state into matlab vector
        //mwSize dims[2];
        //dims[0] = state.rows();
        //dims[1] = 1;
        //mxArray *returnData = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
        //memcpy(mxGetPr(returnData), state.data(), sizeof(double)*state.rows());
        //plhs[0] = returnData;
        return;
    }
    
    //set velocity level DOFs
    if (!strcmp("setQDot", cmd)) {
        
        // Check parameters
        if (nlhs < 0 || nrhs < 3)
            mexErrMsgTxt("Not enough parameters to state.");
        
        double* A = 0;
        size_t  m = 0;
        size_t  n = 0;
        
        A = mxGetPr(prhs[2]);
        // get the dimensions of the first parameter
        m = mxGetM(prhs[2]);
        n = mxGetN(prhs[2]);
        
        Eigen::Map<Eigen::VectorXd> state = mapStateEigen<1>(*dummy_instance);
        
        if(m != state.rows()) {
            mexPrintf("Input vector is of size %i x %i but state is of size %i \n", m,n, state.rows(), 1);
        }
        
        state = Eigen::Map<Eigen::VectorXd>(A, state.rows());
        
        //copy state into matlab vector
        //mwSize dims[2];
        //dims[0] = state.rows();
        //dims[1] = 1;
        //mxArray *returnData = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
        //memcpy(mxGetPr(returnData), state.data(), sizeof(double)*state.rows());
        //plhs[0] = returnData;
        return;
    }
    
    // get mass matrix of the world
    if (!strcmp("M", cmd)) {
        // Check parameters
        if (nlhs < 0 || nrhs < 2)
            mexErrMsgTxt("M: not enough arguments.");
        // Call the method
        
        //run the assembler, copy sparse matrix into matlab sparse matrix and be done with it
        //build mass and stiffness matrices
        AssemblerEigenSparseMatrix<double> mass;
        getMassMatrix(mass, *dummy_instance);
        
        plhs[0] = eigenSparseToMATLAB(*mass);
        return;
    }
    
    // get stiffness matrix of the world
    if (!strcmp("K", cmd)) {
        // Check parameters
        if (nlhs < 0 || nrhs < 2)
            mexErrMsgTxt("K: not enough arguments.");
        // Call the method
        AssemblerEigenSparseMatrix<double> stiffness;
        getStiffnessMatrix(stiffness, *dummy_instance);
        plhs[0] = eigenSparseToMATLAB(*stiffness);
        return;
    }
    
    // get force vector
    if (!strcmp("f", cmd)) {
        // Check parameters
        if (nlhs < 0 || nrhs < 2)
            mexErrMsgTxt("f: not enough arguments.");
        // Call the method
        //dummy_instance->test();
        AssemblerEigenVector<double> force;
        getForceVector(force, *dummy_instance);
        plhs[0] = eigenDenseToMATLAB(*force);
        return;
    }
    
    // get force vector
    if (!strcmp("if", cmd)) {
        // Check parameters
        if (nlhs < 0 || nrhs < 2)
            mexErrMsgTxt("if: not enough arguments.");
        // Call the method
        //dummy_instance->test();
        AssemblerEigenVector<double> force;
        getInternalForceVector(force, *dummy_instance);
        plhs[0] = eigenDenseToMATLAB(*force);
        return;
    }
    // get strain energy
    if (!strcmp("strener", cmd)) {
        // Check parameters
        if (nlhs < 0 || nrhs < 2)
            mexErrMsgTxt("Total Energy: not enough arguments.");
        // Call the method
        plhs[0] = mxCreateNumericMatrix(1,  1, mxDOUBLE_CLASS, mxREAL);
        mxGetPr(plhs[0])[0] = getStrainEnergy(*dummy_instance);
        return;
    }
    
    // get potential done by bodyforces
    if (!strcmp("bdfener", cmd)) {
        // Check parameters
        if (nlhs < 0 || nrhs < 2)
            mexErrMsgTxt("bdfener: not enough arguments.");
        // Call the method
        plhs[0] = mxCreateNumericMatrix(1,  1, mxDOUBLE_CLASS, mxREAL);
        mxGetPr(plhs[0])[0] = getBodyForceEnergy(*dummy_instance);
        return;
    }
    
    
    if (!strcmp("strenertet", cmd)) {
        // Check parameters
        if (nlhs < 0 || nrhs < 2)
            mexErrMsgTxt("strainEnergy: not enough arguments.");
        // Get the per tet strain energy
        Eigen::VectorXd strainEnergy;
        State<double> &state = dummy_instance->getState();
        
        forEach(dummy_instance->getSystemList(), [&strainEnergy, &state](auto a) {
                strainEnergy.resize(impl(a).getF().rows(), 1);
                strainEnergy = impl(a).getStrainEnergyPerElement(state);
        });
        
        plhs[0] = eigenDenseToMATLAB(strainEnergy);

        return;
    }
    
    //get cauchy stresses for all elements
    if(!strcmp("stress",cmd)) {
        if (nlhs < 0 || nrhs < 3)
            mexErrMsgTxt("stress: not enough arguments.");
        // Call the method
        Eigen::MatrixXd stresses;
        Eigen::Matrix<double, 3,3> stress;
        //should find way to avoid this copy
        State<double> state = dummy_instance->getState().mappedState(mxGetPr(prhs[2]));
        
        forEach(dummy_instance->getSystemList(), [&stresses, &stress, &state](auto a) { \
            stresses.resize(impl(a).getF().rows(), 6);
            for(unsigned int ii=0; ii<impl(a).getF().rows(); ++ii) {
                impl(a).getElement(ii)->getCauchyStress(stress, Vec3d(0,0,0), state);
                stresses(ii,0) = stress(0,0);
                stresses(ii,1) = stress(1,1);
                stresses(ii,2) = stress(2,2);
                stresses(ii,3) = stress(1,2);
                stresses(ii,4) = stress(0,2);
                stresses(ii,5) = stress(0,1);
            }
        });

        plhs[0] = eigenDenseToMATLAB(stresses);
        return;
    }
    //need methods for getting mass matrix, stiffness matrix and forces
    
    // Got here, so command not recognized
    mexErrMsgTxt("Command not recognized.");
}
