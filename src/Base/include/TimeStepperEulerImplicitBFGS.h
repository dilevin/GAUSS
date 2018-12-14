//
//  TimeStepperEulerImplicit.h
//  Gauss
//
//  Created by David Levin on 11/23/17.
//
//

#ifndef TimeStepperEulerImplicitBFGS_h
#define TimeStepperEulerImplicitBFGS_h


#include <World.h>
#include <Assembler.h>
#include <TimeStepper.h>
#include <Eigen/Dense>

#include <Eigen/Sparse>
#include <UtilitiesEigen.h>
#include <UtilitiesMATLAB.h>
#include <LBFGS.h>

//TODO Solver Interface
namespace Gauss {
    
    //Given Initial state, step forward in time using linearly implicit Euler Integrator
    template<typename DataType, typename MatrixAssembler, typename VectorAssembler>
    class TimeStepperImplEulerImplicitBFGS
    {
    public:
        
        template<typename Matrix>
        TimeStepperImplEulerImplicitBFGS(Matrix &P, bool precondition=false)  {
            
           //P is constraint projection matrix
            m_P = P;
            m_precondition = precondition;
            
            //setup solver
            LBFGSpp::LBFGSParam<DataType> param;
            param.epsilon = 1e-1;
            param.max_iterations = 1000;
            param.past = 2;
            param.m = 5;
            param.linesearch = LBFGSpp::LBFGS_LINESEARCH_BACKTRACKING_WOLFE;
            
            m_bfgsSolver.setParam(param);
            
        }
        
        TimeStepperImplEulerImplicitBFGS(const TimeStepperImplEulerImplicitBFGS &toCopy) {
            
        }
        
        ~TimeStepperImplEulerImplicitBFGS() { }
        
        //Methods
        //init() //initial conditions will be set at the begining
        template<typename World>
        void step(World &world, double dt, double t);
        
        inline typename VectorAssembler::MatrixType & getLagrangeMultipliers() { return m_lagrangeMultipliers; }
        
    protected:
        
        MatrixAssembler m_massMatrix;
        VectorAssembler m_forceVector;
        MatrixAssembler m_stiffnessMatrix;
        
        Eigen::SparseMatrix<DataType> m_P;
        
        bool m_precondition;
        
        LBFGSpp::LBFGSSolver<DataType> m_bfgsSolver;
        
        typename VectorAssembler::MatrixType m_lagrangeMultipliers;
        
        //#ifdef GAUSS_PARDISO
          //  SolverPardiso<Eigen::SparseMatrix<DataType, Eigen::RowMajor> > m_solver;
        //#else
            Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > m_solver;
        //#endif
        
    private:
    };
}

template<typename DataType, typename MatrixAssembler, typename VectorAssembler>
template<typename World>
void TimeStepperImplEulerImplicitBFGS<DataType, MatrixAssembler, VectorAssembler>::step(World &world, double dt, double t) {
    
    //First two lines work around the fact that C++11 lambda can't directly capture a member variable.
    MatrixAssembler &massMatrix = m_massMatrix;
    MatrixAssembler &stiffnessMatrix = m_stiffnessMatrix;
    VectorAssembler &forceVector = m_forceVector;
    Eigen::SparseMatrix<DataType> &P = m_P;
    
    //Grab the state
    Eigen::VectorXd q = mapStateEigen<0>(world);
    Eigen::VectorXd qDot = mapStateEigen<1>(world);
    
    
    ///get mass matrix
    ASSEMBLEMATINIT(massMatrix, world.getNumQDotDOFs(), world.getNumQDotDOFs());
    ASSEMBLELIST(massMatrix, world.getSystemList(), getMassMatrix);
    ASSEMBLEEND(massMatrix);
    
    //we're going to build equality constraints into our gradient and hessian calcuations
    auto objective = [&world, &q, &qDot, &dt, &forceVector, &massMatrix, &P](Eigen::VectorXd &x, Eigen::VectorXd &grad) -> double {
    
        //std::cout<<"X\n"<<x<<"\n\n";
        mapStateEigen<1>(world) = P.transpose()*x;
        mapStateEigen<0>(world) = q+dt*P.transpose()*x;
        
        ASSEMBLEVECINIT(forceVector, world.getNumQDotDOFs());
        ASSEMBLELIST(forceVector, world.getForceList(), getForce);
        ASSEMBLELIST(forceVector, world.getSystemList(), getForce);
        ASSEMBLEEND(forceVector);
        
        (*forceVector) *= -dt;
        
        (*forceVector) += (*massMatrix)*(P.transpose()*x - qDot);
        
        grad = P*(*forceVector);
        
        double E = (getEnergy(world) - x.transpose()*P*(*massMatrix)*qDot);
    
        return E;
        
    };
    
    double fx = 0.0;
    Eigen::VectorXd qNew = P*qDot;
    qNew.setZero();
   
    if(m_precondition) {
        ASSEMBLEMATINIT(stiffnessMatrix, world.getNumQDotDOFs(), world.getNumQDotDOFs());
        ASSEMBLELIST(stiffnessMatrix, world.getSystemList(), getStiffnessMatrix);
        ASSEMBLEEND(stiffnessMatrix);
        
        Eigen::SparseMatrix<DataType, Eigen::RowMajor> systemMatrix = P*((*m_massMatrix)- dt*dt*(*m_stiffnessMatrix))*P.transpose();
        m_bfgsSolver.minimizeWithPreconditioner(objective, qNew, fx, systemMatrix, m_solver);
    } else {
        m_bfgsSolver.minimize(objective, qNew, fx);
    }
    
    mapStateEigen<1>(world) = P.transpose()*qNew;
    mapStateEigen<0>(world) = q + dt*P.transpose()*qNew;
    
    
}

template<typename DataType, typename MatrixAssembler, typename VectorAssembler>
using TimeStepperEulerImplicitBFGS = TimeStepper<DataType, TimeStepperImplEulerImplicitBFGS<DataType, MatrixAssembler, VectorAssembler> >;


#endif /* TimeStepperEulerImplicit_h */
