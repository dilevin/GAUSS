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
        TimeStepperImplEigenFitSMWImpl(Matrix &P, unsigned int numModes) {
            
            m_numModes = (numModes);
            
            //            std::cout<<m_P.rows()<<std::endl;
            m_P = P;
            m_factored = false;
            m_refactor = true;
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
        
        //        //Ratios diagonal matrix, stored as vector
        Eigen::VectorXd m_R;
        
        
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
        
//        // for calculating the residual. ugly
//        MatrixAssembler m_massMatrixNew;
//        MatrixAssembler m_stiffnessMatrixNew;
//        VectorAssembler m_forceVectorNew;
//        VectorAssembler m_fExtNew;
        
        
        Eigen::SparseMatrix<DataType> m_P;
        
        //storage for lagrange multipliers
        typename VectorAssembler::MatrixType m_lagrangeMultipliers;
        
        bool m_factored, m_refactor;
#ifdef GAUSS_PARDISO
        
        SolverPardiso<Eigen::SparseMatrix<DataType, Eigen::RowMajor> > m_pardiso;
//        SolverPardiso<Eigen::SparseMatrix<DataType, Eigen::RowMajor> > m_pardiso_res;
#else
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
    
//    Eigen::saveMarket(*massMatrix, "massMatrix_woP.dat");

    
    //get stiffness matrix
    ASSEMBLEMATINIT(stiffnessMatrix, world.getNumQDotDOFs(), world.getNumQDotDOFs());
    ASSEMBLELIST(stiffnessMatrix, world.getSystemList(), getStiffnessMatrix);
    ASSEMBLELIST(stiffnessMatrix, world.getForceList(), getStiffnessMatrix);
    ASSEMBLEEND(stiffnessMatrix);
    
//    Eigen::saveMarket(*stiffnessMatrix, "stiffnessMatrix_woP.dat");

    
    //Need to filter internal forces seperately for this applicat
    ASSEMBLEVECINIT(forceVector, world.getNumQDotDOFs());
    ASSEMBLELIST(forceVector, world.getSystemList(), getImpl().getInternalForce);
    ASSEMBLEEND(forceVector);
    
    ASSEMBLEVECINIT(fExt, world.getNumQDotDOFs());
    ASSEMBLELIST(fExt, world.getSystemList(), getImpl().getBodyForce);
    ASSEMBLEEND(fExt);
    
    //constraint Projection
    //    std::cout<<m_P.rows()<<std::endl;
    //    std::cout<<massMatrix.rows()<<std::endl;
    (*massMatrix) = m_P*(*massMatrix)*m_P.transpose();
//    Eigen::saveMarket(*massMatrix, "massMatrix_wP.dat");
    (*stiffnessMatrix) = m_P*(*stiffnessMatrix)*m_P.transpose();
//    Eigen::saveMarket(*stiffnessMatrix, "stiffnessMatrix_wP.dat");

    (*forceVector) = m_P*(*forceVector);
    //Eigendecomposition
    
    // if number of modes not equals to 0, use EigenFitLinear
    if (m_numModes != 0) {
        static_cast<EigenFit*>(std::get<0>(world.getSystemList().getStorage())[0])->calculateEigenFitData(world.getState(),massMatrix,stiffnessMatrix,m_coarseUs,Y,Z);
        
        //    Correct Forces
        (*forceVector) = (*forceVector) + Y*m_coarseUs.first.transpose()*(*forceVector) + m_P*(*fExt);
        
        //Grab the state
        Eigen::Map<Eigen::VectorXd> q = mapStateEigen<0>(world);
        Eigen::Map<Eigen::VectorXd> qDot = mapStateEigen<1>(world);
//        Eigen::VectorXd qDotOld = qDot; // deep copy
        
        
        //setup RHS
        (*forceVector) = (*massMatrix)*m_P*qDot + dt*(*forceVector);
        
        Eigen::VectorXd x0;
        Eigen::SparseMatrix<DataType, Eigen::RowMajor> systemMatrix = (*m_massMatrix)- dt*dt*(*m_stiffnessMatrix);
        
        
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
        
        
        Y = dt*Y;
        Z = -dt*Z;
        
#ifdef GAUSS_PARDISO
        
        m_pardiso.solve(Y);
        Eigen::VectorXd bPrime = Y*(Eigen::MatrixXd::Identity(m_numModes,m_numModes) + Z*m_pardiso.getX()).ldlt().solve(Z*x0);
        
        
#else
        Eigen::VectorXd bPrime = Y*(Eigen::MatrixXd::Identity(m_numModes,m_numModes) + Z*solver.solve(Y)).ldlt().solve(Z*x0);
        
#endif
        
        
#ifdef GAUSS_PARDISO
        
        m_pardiso.solve(bPrime);
        
        x0 -= m_pardiso.getX();
        
        m_pardiso.cleanup();
#else
        
        x0 -= solver.solve(bPrime);
        
#endif
        
        qDot = m_P.transpose()*x0;
        
        //update state
        q = q + dt*qDot;
        
        
//        // calculate the residual. brute force for now. ugly
//        //First two lines work around the fact that C++11 lambda can't directly capture a member variable.
//        MatrixAssembler &massMatrixNew = m_massMatrixNew;
//        MatrixAssembler &stiffnessMatrixNew = m_stiffnessMatrixNew;
//        VectorAssembler &forceVectorNew = m_forceVectorNew;
//        VectorAssembler &fExtNew = m_fExtNew;
//        
//        //get mass matrix
//        ASSEMBLEMATINIT(massMatrixNew, world.getNumQDotDOFs(), world.getNumQDotDOFs());
//        ASSEMBLELIST(massMatrixNew, world.getSystemList(), getMassMatrix);
//        ASSEMBLEEND(massMatrixNew);
//        
//        
//        //get stiffness matrix
//        ASSEMBLEMATINIT(stiffnessMatrixNew, world.getNumQDotDOFs(), world.getNumQDotDOFs());
//        ASSEMBLELIST(stiffnessMatrixNew, world.getSystemList(), getStiffnessMatrix);
//        ASSEMBLELIST(stiffnessMatrixNew, world.getForceList(), getStiffnessMatrix);
//        ASSEMBLEEND(stiffnessMatrixNew);
//        
//        //Need to filter internal forces seperately for this applicat
//        ASSEMBLEVECINIT(forceVectorNew, world.getNumQDotDOFs());
//        ASSEMBLELIST(forceVectorNew, world.getSystemList(), getImpl().getInternalForce);
//        ASSEMBLEEND(forceVectorNew);
//        
//        ASSEMBLEVECINIT(fExtNew, world.getNumQDotDOFs());
//        ASSEMBLELIST(fExtNew, world.getSystemList(), getImpl().getBodyForce);
//        ASSEMBLEEND(fExtNew);
//        
//        //constraint Projection
//        //    std::cout<<m_P.rows()<<std::endl;
//        //    std::cout<<massMatrix.rows()<<std::endl;
//        (*massMatrixNew) = m_P*(*massMatrixNew)*m_P.transpose();
//        (*stiffnessMatrixNew) = m_P*(*stiffnessMatrixNew)*m_P.transpose();
//        (*forceVectorNew) = m_P*(*forceVectorNew);
//        
//        static_cast<EigenFit*>(std::get<0>(world.getSystemList().getStorage())[0])->calculateEigenFitData(world.getState(),massMatrixNew,stiffnessMatrixNew,m_coarseUs,Y,Z);
//        
//        //    Correct Forces
//        (*forceVectorNew) = (*forceVectorNew) + Y*m_coarseUs.first.transpose()*(*forceVectorNew) + m_P*(*fExtNew);
//        
////        //setup RHS
////        (*forceVectorNew) = (*massMatrixNew)*m_P*qDot + dt*(*forceVectorNew);
////
//#ifdef GAUSS_PARDISO
//        // only works for pardiso now
//        m_pardiso_res.symbolicFactorization(*massMatrixNew);
//        m_pardiso_res.numericalFactorization();
//        m_pardiso_res.solve(*forceVectorNew);
//        
//        std::ofstream ofile;
//        ofile.open("residualSI.txt", std::ios::app); //app is append which means it will put the text at the end
//        ofile << (m_P* (qDot - qDotOld) - dt * (m_pardiso_res.getX())).norm() << std::endl;
//        ofile.close();
//#endif
    }
    else // number of modes is 0, so don't need to call EigenFitLinear
    {
        //Grab the state
        Eigen::Map<Eigen::VectorXd> q = mapStateEigen<0>(world);
        Eigen::Map<Eigen::VectorXd> qDot = mapStateEigen<1>(world);
//        Eigen::VectorXd qDotOld = qDot; // deep copy
        
        //    add external Forces
        (*forceVector) = (*forceVector) + m_P*(*fExt);
        
        //setup RHS
        (*forceVector) = (*massMatrix)*m_P*qDot + dt*(*forceVector);
        
        Eigen::VectorXd x0;
        Eigen::SparseMatrix<DataType, Eigen::RowMajor> systemMatrix = (*m_massMatrix)- dt*dt*(*m_stiffnessMatrix);
        
        
        
#ifdef GAUSS_PARDISO
        if(m_refactor || !m_factored) {
            m_pardiso.symbolicFactorization(systemMatrix);
            m_pardiso.numericalFactorization();
            m_factored = true;
        }
        
        m_pardiso.solve(*forceVector);
        x0 = m_pardiso.getX();
        //m_pardiso.cleanup();
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
        
        qDot = m_P.transpose()*x0;
        //    qDot = x0.head(world.getNumQDotDOFs());
        
        //    m_lagrangeMultipliers = x0.tail(world.getNumConstraints());
        
        //update state
        //    updateState(world, world.getState(), dt);
        
        q = q + dt*qDot;
        
//
//        // calculate residual
//        //First two lines work around the fact that C++11 lambda can't directly capture a member variable.
//        //First two lines work around the fact that C++11 lambda can't directly capture a member variable.
//        MatrixAssembler &massMatrixNew = m_massMatrixNew;
//        MatrixAssembler &stiffnessMatrixNew = m_stiffnessMatrixNew;
//        VectorAssembler &forceVectorNew = m_forceVectorNew;
//        VectorAssembler &fExtNew = m_fExtNew;
//
//        //get mass matrix
//        ASSEMBLEMATINIT(massMatrixNew, world.getNumQDotDOFs(), world.getNumQDotDOFs());
//        ASSEMBLELIST(massMatrixNew, world.getSystemList(), getMassMatrix);
//        ASSEMBLEEND(massMatrixNew);
//
//
//        //get stiffness matrix
//        ASSEMBLEMATINIT(stiffnessMatrixNew, world.getNumQDotDOFs(), world.getNumQDotDOFs());
//        ASSEMBLELIST(stiffnessMatrixNew, world.getSystemList(), getStiffnessMatrix);
//        ASSEMBLELIST(stiffnessMatrixNew, world.getForceList(), getStiffnessMatrix);
//        ASSEMBLEEND(stiffnessMatrixNew);
//
//        //Need to filter internal forces seperately for this applicat
//        ASSEMBLEVECINIT(forceVectorNew, world.getNumQDotDOFs());
//        ASSEMBLELIST(forceVectorNew, world.getSystemList(), getImpl().getInternalForce);
//        ASSEMBLEEND(forceVectorNew);
//
//        ASSEMBLEVECINIT(fExtNew, world.getNumQDotDOFs());
//        ASSEMBLELIST(fExtNew, world.getSystemList(), getImpl().getBodyForce);
//        ASSEMBLEEND(fExtNew);
//
//        (*massMatrixNew) = m_P*(*massMatrixNew)*m_P.transpose();
//        (*stiffnessMatrixNew) = m_P*(*stiffnessMatrixNew)*m_P.transpose();
//        (*forceVectorNew) = m_P*(*forceVectorNew);
//
//        //    add external Forces
//        (*forceVectorNew) = (*forceVectorNew) + m_P*(*fExtNew);
//
//#ifdef GAUSS_PARDISO
//        // only works for pardiso now
//        m_pardiso_res.symbolicFactorization(*massMatrixNew);
//        m_pardiso_res.numericalFactorization();
//        m_pardiso_res.solve(*forceVectorNew);
//
//        std::ofstream ofile;
//        ofile.open("residualSI.txt", std::ios::app); //app is append which means it will put the text at the end
//            ofile << (m_P* (qDot - qDotOld) - dt * (m_pardiso_res.getX())).norm() << std::endl;
//        ofile.close();
//#endif
    }
    
    

//    
//#ifdef EDWIN_DEBUG
//    std::ofstream ofile;
//    ofile.open("KE.txt", std::ios::app); //app is append which means it will put the text at the end
//    ofile << std::get<0>(world.getSystemList().getStorage())[0]->getImpl().getKineticEnergy(world.getState()) << std::endl;
//    ofile.close();
////    unsigned int file_ind = 0;
////    std::string name = "pos";
////    std::string fformat = ".obj";
////    std::string filename = name + std::to_string(file_ind) + fformat;
////    struct stat buf;
////    while (stat(filename.c_str(), &buf) != -1)
////    {
////        file_ind++;
////        filename = name + std::to_string(file_ind) + fformat;
////    }
//    unsigned int file_ind = 0;
//    std::string name = "pos";
//    std::string fformat = ".obj";
//    std::string filename = name + std::to_string(file_ind) + fformat;
//    struct stat buf;
//    while (stat(filename.c_str(), &buf) != -1)
//    {
//        file_ind++;
//        filename = name + std::to_string(file_ind) + fformat;
//    }
//
////    unsigned int idx;
////    idx = 0;
////    // getGeometry().first is V
////    Eigen::MatrixXd V_disp = std::get<0>(world.getSystemList().getStorage())[0]->getGeometry().first;
////    
////    for(unsigned int vertexId=0;  vertexId < std::get<0>(world.getSystemList().getStorage())[0]->getGeometry().first.rows(); ++vertexId) {
////        
////        // because getFinePosition is in EigenFit, not another physical system Impl, so don't need getImpl()
////        V_disp(vertexId,0) += q(idx);
////        idx++;
////        V_disp(vertexId,1) += q(idx);
////        idx++;
////        V_disp(vertexId,2) += q(idx);
////        idx++;
////    }
////    igl::writeOBJ(filename,V_disp,std::get<0>(world.getSystemList().getStorage())[0]->getGeometry().second);
////    saveMarket(V_disp, "V.dat");
////    saveMarket(std::get<0>(world.getSystemList().getStorage())[0]->getGeometry().second, "F.dat");
//#endif
}

template<typename DataType, typename MatrixAssembler, typename VectorAssembler>
using TimeStepperEigenFitSMW = TimeStepper<DataType, TimeStepperImplEigenFitSMWImpl<DataType, MatrixAssembler, VectorAssembler> >;




#endif /* TimeStepperEigenFitSMW_h */
