#include "mex.h"
#include "class_handle.hpp"

//Gauss Stuff
#include <GaussIncludes.h>
#include <FEMIncludes.h>
#include <ForceSpring.h>

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

//Generic world object that supports the sim objects I need
typedef World<  double,
                std::tuple<
                    FEMLinearTets *,
                    FEMPlaneStrainTri *
                >,
                std::tuple<ForceSpring<double> *>,
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
        
        char * femType = mxArrayToString(prhs[0]);
        
        if(strcmp(femType, "elastic_linear_tetrahedra") == 0) {
            mexPrintf("Initialize Linear Elastic Tetrahedra");
            FEMLinearTets *FEM = new FEMLinearTets(matlabToDouble(prhs[2]), matlabToInt32(prhs[3]));
            world->addSystem(FEM);
        } else if(femType, "elastic_linear_plane_strain_tri") {
            mexPrintf("Initialize Linear Elastic Plane Strain Tri\n");
            FEMPlaneStrainTri *FEM = new FEMPlaneStrainTri(matlabToDouble(prhs[2]), matlabToInt32(prhs[3]));
            world->addSystem(FEM);
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
        if (nlhs < 1 || nrhs < 1)
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
    
    //need methods for getting mass matrix, stiffness matrix and forces
    
    // Got here, so command not recognized
    mexErrMsgTxt("Command not recognized.");
}
