//
//  TimeStepperEulerImplicitLinearCProjection.h
//  Gauss
//
//  Created by Yu Ju Edwin Chen on 2018-06-18.
//

#ifndef TimeStepperEulerImplicitLinearCProjection_h
#define TimeStepperEulerImplicitLinearCProjection_h


#include <World.h>
#include <Assembler.h>
#include <TimeStepper.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <UtilitiesEigen.h>
#include <UtilitiesMATLAB.h>
#include <Eigen/SparseCholesky>
#include <SolverPardiso.h>

// use constraint projection instead of using KKT. similar to TimeStepperEigenFitSMW

//TODO Solver Interface
namespace Gauss {
    
    //Given Initial state, step forward in time using linearly implicit Euler Integrator
    template<typename DataType, typename MatrixAssembler, typename VectorAssembler>
    class TimeStepperImplEulerImplicitLinearCProjection
    {
    public:
        
//        TimeStepperImplEulerImplicitLinear(bool refactor = true) {
//            m_factored = false;
//            m_refactor = refactor;
//        }
//
//        TimeStepperImplEulerImplicitLinear(const TimeStepperImplEulerImplicitLinear &toCopy) {
//            m_factored = false;
//        }
        
        template<typename Matrix>
        TimeStepperImplEulerImplicitLinearCProjection(Matrix &P) {
            
            m_P = P;
            m_factored = false;
            m_refactor = false;
        }

        
        ~TimeStepperImplEulerImplicitLinearCProjection() {
#ifdef GAUSS_PARDISO
            m_pardiso.cleanup();
#endif
        }
        
        //Methods
        template<typename World>
        void step(World &world, double dt, double t);
        
        inline typename VectorAssembler::MatrixType & getLagrangeMultipliers() { return m_lagrangeMultipliers; }
        
    protected:
        
        MatrixAssembler m_massMatrix;
        MatrixAssembler m_stiffnessMatrix;
        VectorAssembler m_forceVector;
        
        Eigen::SparseMatrix<DataType> m_P;
        
        //storage for lagrange multipliers
        typename VectorAssembler::MatrixType m_lagrangeMultipliers;
        
        bool m_factored, m_refactor;
        
#ifdef GAUSS_PARDISO
        
        SolverPardiso<Eigen::SparseMatrix<DataType, Eigen::RowMajor> > m_pardiso;
        
#endif
        
    private:
    };
}

template<typename DataType, typename MatrixAssembler, typename VectorAssembler>
template<typename World>
void TimeStepperImplEulerImplicitLinearCProjection<DataType, MatrixAssembler, VectorAssembler>::step(World &world, double dt, double t) {
    
    
    Eigen::SparseMatrix<DataType, Eigen::RowMajor> systemMatrix;
    Eigen::VectorXd x0;
    
//    if(m_refactor || !m_factored) {
    
        //First two lines work around the fact that C++11 lambda can't directly capture a member variable.
        MatrixAssembler &massMatrix = m_massMatrix;
        MatrixAssembler &stiffnessMatrix = m_stiffnessMatrix;
        
        //get mass matrix
        ASSEMBLEMATINIT(massMatrix, world.getNumQDotDOFs(), world.getNumQDotDOFs());
        ASSEMBLELIST(massMatrix, world.getSystemList(), getMassMatrix);
        
//        //add in the constraints
//        ASSEMBLELISTOFFSET(massMatrix, world.getConstraintList(), getGradient, world.getNumQDotDOFs(), 0);
//        ASSEMBLELISTOFFSETTRANSPOSE(massMatrix, world.getConstraintList(), getGradient, 0, world.getNumQDotDOFs());
//        ASSEMBLEEND(massMatrix);
    
        
        //get stiffness matrix
        ASSEMBLEMATINIT(stiffnessMatrix, world.getNumQDotDOFs(), world.getNumQDotDOFs());
        ASSEMBLELIST(stiffnessMatrix, world.getSystemList(), getStiffnessMatrix);
        ASSEMBLELIST(stiffnessMatrix, world.getForceList(), getStiffnessMatrix);
        ASSEMBLEEND(stiffnessMatrix);
    
    // constraint projection only for fixed or jump
    (*massMatrix) = m_P*(*massMatrix)*m_P.transpose();
    (*stiffnessMatrix) = m_P*(*stiffnessMatrix)*m_P.transpose();
    

        systemMatrix = (*m_massMatrix)- dt*dt*(*m_stiffnessMatrix);
        
//    }
    
    VectorAssembler &forceVector = m_forceVector;
    
    ASSEMBLEVECINIT(forceVector, world.getNumQDotDOFs());
    ASSEMBLELIST(forceVector, world.getForceList(), getForce);
    ASSEMBLELIST(forceVector, world.getSystemList(), getForce);
//    ASSEMBLELISTOFFSET(forceVector, world.getConstraintList(), getDbDt, world.getNumQDotDOFs(), 0);
    ASSEMBLEEND(forceVector);
    
    // constraint projection only for fixed or jump
    (*forceVector) = m_P*(*forceVector);
    
    //Grab the state
    Eigen::Map<Eigen::VectorXd> q = mapStateEigen<0>(world);
    Eigen::Map<Eigen::VectorXd> qDot = mapStateEigen<1>(world);
    
    //setup RHS
//    (*forceVector).head(world.getNumQDotDOFs()) = (*m_massMatrix).block(0,0, world.getNumQDotDOFs(), world.getNumQDotDOFs())*qDot + dt*(*forceVector).head(world.getNumQDotDOFs());
    (*forceVector) = (*massMatrix)*m_P*qDot + dt*(*forceVector);

    
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
}

template<typename DataType, typename MatrixAssembler, typename VectorAssembler>
using TimeStepperEulerImplicitLinearCProjection = TimeStepper<DataType, TimeStepperImplEulerImplicitLinearCProjection<DataType, MatrixAssembler, VectorAssembler> >;





#endif /* TimeStepperEulerImplicitLinearCProjection_h */
