//
//  TimeStepperEigenFitSMW.h
//  Gauss
//
//  Created by Edwin Chen on 2018-05-15.
//
//

#ifndef TimeStepperEigenFitSMW_h
#define TimeStepperEigenFitSMW_h



#include <World.h>
#include <Assembler.h>
#include <TimeStepper.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <UtilitiesEigen.h>
#include <UtilitiesMATLAB.h>
#include <Eigen/SparseCholesky>
#include <SolverPardiso.h>
#include <EigenFit.h>

//TODO Solver Interface
namespace Gauss {
    
    //Given Initial state, step forward in time using linearly implicit Euler Integrator
    template<typename DataType, typename MatrixAssembler, typename VectorAssembler>
    class TimeStepperImplEigenFitSMWImpl
    {
    public:
        
        template<typename Matrix>
        TimeStepperImplEigenFitSMWImpl(Matrix &P) {
            
            m_numToCorrect = 10;
            m_numModes = 10;
//            m_R.setConstant(m_numModes, 1.0);
//            m_I.setConstant(m_numModes, 1.0);
            m_P = P;
        }
        
        TimeStepperImplEigenFitSMWImpl(const TimeStepperImplEigenFitSMWImpl &toCopy) {
            
        }
        
        ~TimeStepperImplEigenFitSMWImpl() { }
        
        //Methods
        //init() //initial conditions will be set at the begining
        template<typename World>
        void step(World &world, double dt, double t);
        
        inline typename VectorAssembler::MatrixType & getLagrangeMultipliers() { return m_lagrangeMultipliers; }
        
    protected:
        
        //num modes to correct
        unsigned int m_numModes;
        unsigned int m_numToCorrect;
        
//        //Ratios diagonal matrix, stored as vector
        Eigen::VectorXd m_R;
//        Eigen::VectorXd m_I;
        
        
        //Subspace Eigenvectors and eigenvalues from the embedded fine mesh
        std::pair<Eigen::MatrixXx<DataType>, Eigen::VectorXx<DataType> > m_fineUs;
        //Subspace Eigenvectors and eigenvalues from this coarse mesh
        std::pair<Eigen::MatrixXx<DataType>, Eigen::VectorXx<DataType> > m_coarseUs;
        // for SMW
        Eigen::MatrixXd Y;
        Eigen::MatrixXd Z;

        MatrixAssembler m_massMatrix;
        MatrixAssembler m_stiffnessMatrix;
        VectorAssembler m_forceVector;
        VectorAssembler m_fExt;
        
        Eigen::SparseMatrix<DataType> m_P;
        
        //storage for lagrange multipliers
        typename VectorAssembler::MatrixType m_lagrangeMultipliers;
        
#ifdef GAUSS_PARDISO
        
        SolverPardiso<Eigen::SparseMatrix<DataType, Eigen::RowMajor> > m_pardiso;
#else
        bool m_factored, m_refactor;
        
//        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > &solver = m_eigensolver;
//        Eigen::VectorXd eigenSolverResult;
#endif
        
    private:
    };
}

template<typename DataType, typename MatrixAssembler, typename VectorAssembler>
template<typename World>
void TimeStepperImplEigenFitSMWImpl<DataType, MatrixAssembler, VectorAssembler>::step(World &world, double dt, double t) {
    
    //First two lines work around the fact that C++11 lambda can't directly capture a member variable.
    MatrixAssembler &massMatrix = m_massMatrix;
    MatrixAssembler &stiffnessMatrix = m_stiffnessMatrix;
    VectorAssembler &forceVector = m_forceVector;
    VectorAssembler &fExt = m_fExt;
    
    //get mass matrix
    ASSEMBLEMATINIT(massMatrix, world.getNumQDotDOFs(), world.getNumQDotDOFs());
    ASSEMBLELIST(massMatrix, world.getSystemList(), getMassMatrix);
    ASSEMBLEEND(massMatrix);
    
    
    //get stiffness matrix
    ASSEMBLEMATINIT(stiffnessMatrix, world.getNumQDotDOFs(), world.getNumQDotDOFs());
    ASSEMBLELIST(stiffnessMatrix, world.getSystemList(), getStiffnessMatrix);
    ASSEMBLELIST(stiffnessMatrix, world.getForceList(), getStiffnessMatrix);
    ASSEMBLEEND(stiffnessMatrix);
    
    //Need to filter internal forces seperately for this applicat
    ASSEMBLEVECINIT(forceVector, world.getNumQDotDOFs());
    ASSEMBLELIST(forceVector, world.getSystemList(), getImpl().getInternalForce);
    ASSEMBLEEND(forceVector);
    
    ASSEMBLEVECINIT(fExt, world.getNumQDotDOFs());
    ASSEMBLELIST(fExt, world.getSystemList(), getImpl().getBodyForce);
    ASSEMBLEEND(fExt);
    
    //constraint Projection
    (*massMatrix) = m_P*(*massMatrix)*m_P.transpose();
    (*stiffnessMatrix) = m_P*(*stiffnessMatrix)*m_P.transpose();
    (*forceVector) = m_P*(*forceVector);
    
//    saveMarket(*massMatrix, "coarsemass.mtx");
    
    //Eigendecomposition
//    m_coarseUs = generalizedEigenvalueProblem((*stiffnessMatrix), (*massMatrix), m_numModes, 0.0);
  
    static_cast<EigenFit*>(std::get<0>(world.getSystemList().getStorage())[0])->calculateEigenFitData(world.getState(),massMatrix,stiffnessMatrix,m_coarseUs,Y,Z);
    
//    Correct Forces
    (*forceVector) = (*forceVector) + Y*m_coarseUs.first.transpose()*(*forceVector) + m_P*(*fExt);
    
    //Grab the state
    Eigen::Map<Eigen::VectorXd> q = mapStateEigen<0>(world);
    Eigen::Map<Eigen::VectorXd> qDot = mapStateEigen<1>(world);
    
    
    //std::cout<<"F: "<<(*forceVector)<<"\n";
    //setup RHS
    (*forceVector) = (*massMatrix)*m_P*qDot + dt*(*forceVector);
    
    Eigen::VectorXd x0;
    Eigen::SparseMatrix<DataType, Eigen::RowMajor> systemMatrix = (*m_massMatrix)- dt*dt*(*m_stiffnessMatrix);
    
//    int file_ind = 0;
//    std::string name = "coarsestiffness_unmodified";
//    std::string fformat = ".dat";
//    std::string filename = name + std::to_string(file_ind) + fformat;
//    struct stat buf;
//    while (stat(filename.c_str(), &buf) != -1)
//    {
//        file_ind++;
//        filename = name + std::to_string(file_ind) + fformat;
//    }
//    saveMarket(*m_stiffnessMatrix, filename);

#ifdef GAUSS_PARDISO

    m_pardiso.symbolicFactorization(systemMatrix, m_numModes);
    m_pardiso.numericalFactorization();
    
//    SMW update for Eigenfit here

    m_pardiso.solve(*forceVector);
    x0 = m_pardiso.getX();

#else
    //solve system (Need interface for solvers but for now just use Eigen LLt)
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
    
    if(m_refactor || !m_factored) {
        solver.compute(systemMatrix);
    }
    
    if(solver.info()!=Eigen::Success) {
        // decomposition failed
        assert(1 == 0);
        std::cout<<"Decomposition Failed \n";
        exit(1);
    }
    
    if(solver.info()!=Eigen::Success) {
        // solving failed
        assert(1 == 0);
        std::cout<<"Solve Failed \n";
        exit(1);
    }
    
    
    x0 = solver.solve((*forceVector));

#endif

//    file_ind = 0;
//    name = "x0_original";
//    fformat = ".dat";
//    filename = name + std::to_string(file_ind) + fformat;
//    while (stat(filename.c_str(), &buf) != -1)
//    {
//        file_ind++;
//        filename = name + std::to_string(file_ind) + fformat;
//    }
//    saveMarket(x0, filename);

    
    Y = dt*Y;
    Z = -dt*Z;
    
#ifdef GAUSS_PARDISO

    m_pardiso.solve(Y);
    Eigen::VectorXd bPrime = Y*(Eigen::MatrixXd::Identity(m_numModes,m_numModes) + Z*m_pardiso.getX()).ldlt().solve(Z*x0);
    
    
#else
    Eigen::VectorXd bPrime = Y*(Eigen::MatrixXd::Identity(m_numModes,m_numModes) + Z*solver.solve(Y)).ldlt().solve(Z*x0);
    
#endif

//    file_ind = 0;
//    name = "bPrime";
//    fformat = ".dat";
//    filename = name + std::to_string(file_ind) + fformat;
//
//    while (stat(filename.c_str(), &buf) != -1)
//    {
//        file_ind++;
//        filename = name + std::to_string(file_ind) + fformat;
//    }
//    saveMarket(bPrime, filename);
//
    
#ifdef GAUSS_PARDISO

    m_pardiso.solve(bPrime);
    
    x0 -= m_pardiso.getX();
    
    m_pardiso.cleanup();
#else
    
    x0 -= solver.solve(bPrime);
    
#endif
//    file_ind = 0;
//    name = "x0_SMW";
//    fformat = ".dat";
//    filename = name + std::to_string(file_ind) + fformat;
//    while (stat(filename.c_str(), &buf) != -1)
//    {
//        file_ind++;
//        filename = name + std::to_string(file_ind) + fformat;
//    }
//    saveMarket(x0, filename);
//
    
    qDot = m_P.transpose()*x0;
    
    //update state
    q = q + dt*qDot;
    
//    file_ind = 0;
//    name = "coarsemesh_q";
//    fformat = ".dat";
//    filename = name + std::to_string(file_ind) + fformat;
//    while (stat(filename.c_str(), &buf) != -1)
//    {
//        file_ind++;
//        filename = name + std::to_string(file_ind) + fformat;
//    }
//    saveMarket(q, filename);

    //std::cout<<"Q: "<<q<<"\n";
}

template<typename DataType, typename MatrixAssembler, typename VectorAssembler>
using TimeStepperEigenFitSMW = TimeStepper<DataType, TimeStepperImplEigenFitSMWImpl<DataType, MatrixAssembler, VectorAssembler> >;




#endif /* TimeStepperEigenFitSMW_h */
