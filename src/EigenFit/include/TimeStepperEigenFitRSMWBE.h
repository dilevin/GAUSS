//
//  TimeStepperEigenFitRSMWBE.h
//  Gauss
//
//  Created by Edwin Chen on 2018-05-15.
//
//

#ifndef TimeStepperEigenFitRSMWBE_h
#define TimeStepperEigenFitRSMWBE_h



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
#include <limits>


//TODO Solver Interface
namespace Gauss {
    
    //Given Initial state, step forward in time using linearly implicit Euler Integrator
    template<typename DataType, typename MatrixAssembler, typename VectorAssembler>
    class TimeStepperImplEigenFitRSMWBEImpl
    {
    public:
        
        template<typename Matrix>
        TimeStepperImplEigenFitRSMWBEImpl(Matrix &P, unsigned int numModes) {
            
            m_numModes = (numModes);
            
            //            std::cout<<m_P.rows()<<std::endl;
            m_P = P;
            m_factored = false;
            // refactor for every solve
            m_refactor = true;
            
            // init residual
            res = std::numeric_limits<double>::infinity();
            
            it_outer = 0;
            it_inner = 0;
            // constants from Nocedal and Wright
            step_size = 1;
            c1 = 1e-4;
            c2 = 0.9;
            
            this->a = 0.0;
            this->b = -0.02;
            
            //            // construct the selection matrix
            int m = P.cols() - P.rows();
            int n = P.cols();
            Eigen::VectorXd one = Eigen::VectorXd::Ones(P.rows());
            Eigen::MatrixXd P_col_sum =   one.transpose() * P;
            std::vector<Eigen::Triplet<DataType> > triplets;
            Eigen::SparseMatrix<DataType> S;
            //            std::cout<<"P_col_sum size: "<< P_col_sum.rows()<< " " << P_col_sum.cols()<<std::endl;
            S.resize(m,n);
            //
            
            unsigned int row_index =0;
            for(unsigned int col_index = 0; col_index < n; col_index++) {
                
                //add triplet into matrix
                //                std::cout<<col_index<<std::endl;
                if (P_col_sum(0,col_index) == 0) {
                    triplets.push_back(Eigen::Triplet<DataType>(row_index,col_index,1));
                    row_index+=1;
                }
                
            }
            
            S.setFromTriplets(triplets.begin(), triplets.end());
            
            m_S = S;
            std::cout<<"m_S sum: "<< m_S.sum()<<std::endl;
            
        }
        
        
        TimeStepperImplEigenFitRSMWBEImpl(const TimeStepperImplEigenFitRSMWBEImpl &toCopy) {
            
        }
        
        ~TimeStepperImplEigenFitRSMWBEImpl() { }
        
        //Methods
        //init() //initial conditions will be set at the begining
        template<typename World>
        void step(World &world, double dt, double t);
        
        inline typename VectorAssembler::MatrixType & getLagrangeMultipliers() { return m_lagrangeMultipliers; }
        
        double* rhs  = NULL;
        double* v_old  = NULL;
        double* v_temp = NULL;
        double* q_old = NULL;
        double* q_temp = NULL;
        
        
        // for damping
        double a;
        // negative for opposite sign for stiffness
        double b;
        
        //        int flag = 0;
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
        Eigen::SparseMatrix<DataType> m_S;
        
        
        //storage for lagrange multipliers
        typename VectorAssembler::MatrixType m_lagrangeMultipliers;
        
        bool m_factored, m_refactor;
        
        // iteration counter
        int it_outer, it_inner;
        // residual
        double res, res_old, step_size, c1, c2;
        
        //        Eigen::VectorXd res;
        
#ifdef GAUSS_PARDISO
        
        SolverPardiso<Eigen::SparseMatrix<DataType, Eigen::RowMajor> > m_pardiso;
        SolverPardiso<Eigen::SparseMatrix<DataType, Eigen::RowMajor> > m_pardiso_test;
        SolverPardiso<Eigen::SparseMatrix<DataType, Eigen::RowMajor> > m_pardiso_mass;
        //        SolverPardiso<Eigen::SparseMatrix<DataType, Eigen::RowMajor> > m_pardiso_res;
#else
#endif
        
    private:
    };
}

template<typename DataType, typename MatrixAssembler, typename VectorAssembler>
template<typename World>
void TimeStepperImplEigenFitRSMWBEImpl<DataType, MatrixAssembler, VectorAssembler>::step(World &world, double dt, double t) {
    
    //First two lines work around the fact that C++11 lambda can't directly capture a member variable.
    MatrixAssembler &massMatrix = m_massMatrix;
    MatrixAssembler &stiffnessMatrix = m_stiffnessMatrix;
    VectorAssembler &forceVector = m_forceVector;
    VectorAssembler &fExt = m_fExt;
    
    // init rhs
    rhs = new double[m_P.rows()];
    Eigen::Map<Eigen::VectorXd> eigen_rhs(rhs,m_P.rows());
    
    // velocity from previous step (for calculating residual)
    v_old = new double[world.getNumQDotDOFs()];
    Eigen::Map<Eigen::VectorXd> eigen_v_old(v_old,world.getNumQDotDOFs());
    
    // velocity from previous step (for calculating residual)
    v_temp = new double[world.getNumQDotDOFs()];
    Eigen::Map<Eigen::VectorXd> eigen_v_temp(v_temp,world.getNumQDotDOFs());
    
    q_old = new double[world.getNumQDotDOFs()];
    Eigen::Map<Eigen::VectorXd> eigen_q_old(q_old,world.getNumQDotDOFs());
    
    q_temp = new double[world.getNumQDotDOFs()];
    Eigen::Map<Eigen::VectorXd> eigen_q_temp(q_temp,world.getNumQDotDOFs());
    
    
    //Grab the state
    Eigen::Map<Eigen::VectorXd> q = mapStateEigen<0>(world);
    Eigen::Map<Eigen::VectorXd> qDot = mapStateEigen<1>(world);
    
    // this is the velocity and q from the end of the last time step
    for (int ind = 0; ind < qDot.rows(); ind++) {
        eigen_v_old(ind) = qDot(ind);
        eigen_q_old(ind) = q(ind);
        eigen_v_temp(ind) = qDot(ind);
    }
    
    
    do {
        std::cout<<"it outer: " << it_outer<<std::endl;
        it_outer = it_outer + 1;
        if (it_outer > 20) {
            std::cout<< "warning: quasi-newton more than 20 iterations." << std::endl;
        }
        
        eigen_q_temp = eigen_q_old + dt * (qDot);
        
        // set the state
        q = eigen_q_temp;
        
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
        
        // Eigen::saveMarket(m_P, "m_P.dat");
        //constraint Projection
        (*massMatrix) = m_P*(*massMatrix)*m_P.transpose();
        (*stiffnessMatrix) = m_P*(*stiffnessMatrix)*m_P.transpose();
        
        (*forceVector) = m_P*(*forceVector);
        //Eigendecomposition
        
        // if number of modes not equals to 0, use EigenFit
        if (m_numModes != 0) {
            
            try{
                if(static_cast<EigenFit*>(std::get<0>(world.getSystemList().getStorage())[0])->calculateEigenFitData(q,massMatrix,stiffnessMatrix,m_coarseUs,Y,Z)) throw 1;
            }catch(...)
            {
                std::cout<<"hausdorff distance check fail\n";
                static_cast<EigenFit*>(std::get<0>(world.getSystemList().getStorage())[0])->flag = 2;
                return;
            }
            
            //    Correct Forces
            (*forceVector) = (*forceVector) + Y*m_coarseUs.first.transpose()*(*forceVector);
            
            
            // add damping
            //            (*forceVector) = (*forceVector) -  (a * (*massMatrix) + b*(*stiffnessMatrix)) * m_P * ( qDot);
            (*forceVector) = (*forceVector) -  (a * (*massMatrix) + b*(*stiffnessMatrix)) * m_P * ( qDot) + b * (Y*(Z * (m_P *( qDot))));
        }
        else
        {
            // add damping
            (*forceVector) = (*forceVector) -  (a * (*massMatrix) + b*(*stiffnessMatrix)) * m_P *(qDot);
        }
        
        // add external force
        (*forceVector) = (*forceVector) + m_P*(*fExt);
        
        
        //setup RHS
        eigen_rhs = (*massMatrix)*m_P*(qDot-eigen_v_old) - dt*(*forceVector);
        
        m_pardiso_mass.symbolicFactorization(*massMatrix);
        m_pardiso_mass.numericalFactorization();
        m_pardiso_mass.solve(*forceVector);
        
        res_old = 1.0/2.0 * dt * dt * ((m_pardiso_mass.getX()).transpose()).squaredNorm();
        
        // std::cout<<"res_old: "<<res_old << std::endl;
        
        Eigen::VectorXd x0;
        // last term is damping
        Eigen::SparseMatrix<DataType, Eigen::RowMajor> systemMatrix = -(*m_massMatrix)+ dt*dt*(*m_stiffnessMatrix)- dt * (a *(*m_massMatrix) + b * (*m_stiffnessMatrix));
        
        
#ifdef GAUSS_PARDISO
        
        m_pardiso.symbolicFactorization(systemMatrix, m_numModes);
        m_pardiso.numericalFactorization();
        
        //    SMW update for Eigenfit here
        
        m_pardiso.solve(eigen_rhs);
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
        
        
        x0 = solver.solve((eigen_rhs));
        
#endif
        
        if (m_numModes != 0) {
            // Y = dt*Y;
            Y = (-(dt*dt) + dt * b)*Y;
            // Y = (dt*dt - dt*b)*Y;
            //
            // Z = dt*Z;
            // todo: add damping here
            
#ifdef GAUSS_PARDISO
            
            m_pardiso.solve(Y);
            Eigen::MatrixXd APrime = Z*m_pardiso.getX();
            Eigen::VectorXd bPrime = Y*(Eigen::MatrixXd::Identity(m_numModes,m_numModes) + APrime).ldlt().solve(Z*x0);
            
            
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
        }
        //        qDot = m_P.transpose()*x0;
        
        auto Dv = m_P.transpose()*x0;
        //        std::cout<<"m_S*Dv"<< m_S*Dv <<std::endl;
        eigen_v_temp = qDot;
        eigen_v_temp = eigen_v_temp + Dv*step_size;
        //update state
        //        std::cout<<"m_S*eigen_v_temp"<< m_S*eigen_v_temp <<std::endl;
        
        q = eigen_q_old + dt*eigen_v_temp;
        
        
        //        std::cout<<"q "<<q.rows()<< std::endl;
        
        // calculate the residual. brute force for now. ugly
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
        //        (*massMatrix) = m_P*(*massMatrix)*m_P.transpose();
        (*stiffnessMatrix) = m_P*(*stiffnessMatrix)*m_P.transpose();
        
        (*forceVector) = m_P*(*forceVector);
        
        if (m_numModes != 0) {
            
            static_cast<EigenFit*>(std::get<0>(world.getSystemList().getStorage())[0])->calculateEigenFitData(q,massMatrix,stiffnessMatrix,m_coarseUs,Y,Z);
            
            //    Correct Forces
            (*forceVector) = (*forceVector) + Y*m_coarseUs.first.transpose()*(*forceVector);
            // add damping
//            (*forceVector) = (*forceVector) -  (a * (*massMatrix) + b*(*stiffnessMatrix)) * m_P * (qDot);
             (*forceVector) = (*forceVector) -  (a * (*massMatrix) + b*(*stiffnessMatrix)) * m_P * (qDot) - b * (Y * (Z * (m_P *(qDot))));
        }
        else
        {
            (*forceVector) = (*forceVector) -  (a * (*massMatrix) + b*(*stiffnessMatrix)) * m_P *(qDot);
            
        }
        // add external force
        (*forceVector) = (*forceVector) + m_P*(*fExt);
        
        m_pardiso_mass.solve(*forceVector);
        std::cout << "res: " << 1.0/2.0 * (m_P*(eigen_v_temp - eigen_v_old) - dt*m_pardiso_mass.getX()).squaredNorm()<< std::endl;
        res  = 1.0/2.0 * (m_P*(eigen_v_temp - eigen_v_old) - dt*m_pardiso_mass.getX()).squaredNorm();
        
        //            while (res > res_old + c1 * step_size * res * 2) {
        //                res_old = res;
        //                it_inner = it_inner + 1;
        //                if (it_inner > 40) {
        //                    std::cout<<"warning: line search more than 40 iterations." << std::endl;
        //                }
        //                step_size = 0.9 * step_size;
        //                eigen_v_temp = qDot + step_size*Dv;
        //                q = eigen_q_old + dt*eigen_v_temp;
        //
        //                // calculate the residual. brute force for now. ugly
        //                //get stiffness matrix
        //                ASSEMBLEMATINIT(stiffnessMatrix, world.getNumQDotDOFs(), world.getNumQDotDOFs());
        //                ASSEMBLELIST(stiffnessMatrix, world.getSystemList(), getStiffnessMatrix);
        //                ASSEMBLELIST(stiffnessMatrix, world.getForceList(), getStiffnessMatrix);
        //                ASSEMBLEEND(stiffnessMatrix);
        //
        //                //Need to filter internal forces seperately for this applicat
        //                ASSEMBLEVECINIT(forceVector, world.getNumQDotDOFs());
        //                ASSEMBLELIST(forceVector, world.getSystemList(), getImpl().getInternalForce);
        //                ASSEMBLEEND(forceVector);
        //
        //                ASSEMBLEVECINIT(fExt, world.getNumQDotDOFs());
        //                ASSEMBLELIST(fExt, world.getSystemList(), getImpl().getBodyForce);
        //                ASSEMBLEEND(fExt);
        //
        //                //constraint Projection
        //                //        (*massMatrix) = m_P*(*massMatrix)*m_P.transpose();
        //                (*stiffnessMatrix) = m_P*(*stiffnessMatrix)*m_P.transpose();
        //
        //                (*forceVector) = m_P*(*forceVector);
        //                static_cast<EigenFit*>(std::get<0>(world.getSystemList().getStorage())[0])->calculateEigenFitData(q,massMatrix,stiffnessMatrix,m_coarseUs,Y,Z);
        //
        //                //    Correct Forces
        //                (*forceVector) = (*forceVector) + Y*m_coarseUs.first.transpose()*(*forceVector);
        //
        //                // add external force
        //                (*forceVector) = (*forceVector) + m_P*(*fExt);
        //
        //                m_pardiso_mass.solve(*forceVector);
        //                std::cout << "res: " << 1.0/2.0 * (m_P*(eigen_v_temp - eigen_v_old) - dt*m_pardiso_mass.getX()).squaredNorm()<< std::endl;
        //                res  = 1.0/2.0 * (m_P*(eigen_v_temp - eigen_v_old) - dt*m_pardiso_mass.getX()).squaredNorm();
        //
        //            }
        
        qDot = eigen_v_temp;
        q = eigen_q_old + dt * (qDot);
        
    } while(res > 1e-6);
    it_outer = 0;
    //    std::cout<<"m_S*qDot"<< m_S*qDot <<std::endl;
    
    q = eigen_q_old + dt*qDot;
}

template<typename DataType, typename MatrixAssembler, typename VectorAssembler>
using TimeStepperEigenFitRSMWBE = TimeStepper<DataType, TimeStepperImplEigenFitRSMWBEImpl<DataType, MatrixAssembler, VectorAssembler> >;




#endif /* TimeStepperEigenFitRSMWBE_h */
