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
        TimeStepperImplEulerImplicitBFGS(Matrix &P)  {
            
           //P is constraint projection matrix
            m_P = P;
            
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
        Eigen::SparseMatrix<DataType> m_P;
        
        typename VectorAssembler::MatrixType m_lagrangeMultipliers;
        
    private:
    };
}

template<typename DataType, typename MatrixAssembler, typename VectorAssembler>
template<typename World>
void TimeStepperImplEulerImplicitBFGS<DataType, MatrixAssembler, VectorAssembler>::step(World &world, double dt, double t) {
    
    //First two lines work around the fact that C++11 lambda can't directly capture a member variable.
    MatrixAssembler &massMatrix = m_massMatrix;
    VectorAssembler &forceVector = m_forceVector;
    Eigen::SparseMatrix<DataType> &P = m_P;
    
    //Grab the state
    Eigen::VectorXd q = mapStateEigen<0>(world);
    Eigen::VectorXd qDot = mapStateEigen<1>(world);
    
    
    ///get mass matrix
    ASSEMBLEMATINIT(massMatrix, world.getNumQDotDOFs(), world.getNumQDotDOFs());
    ASSEMBLELIST(massMatrix, world.getSystemList(), getMassMatrix);
    ASSEMBLEEND(massMatrix);
    
    //std::cout<<"F: "<<(*forceVector)<<"\n";
    //setup RHS
    
    
    //toMatlab(*massMatrix, "./testMHex.txt");
    //toMatlab(*stiffnessMatrix, "./testKHex.txt");
    //toMatlab(*forceVector, "./testfHex.txt");
    
    //std::cout<<"F: \n"<<(*forceVector)<<"\n";
    
    
    //Eigen::SparseMatrix<DataType, Eigen::RowMajor> systemMatrix = (*m_massMatrix)- dt*dt*(*m_stiffnessMatrix);
    
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
    
    LBFGSpp::LBFGSParam<DataType> param;
    param.epsilon = 1e-1;
    param.max_iterations = 1000;
    param.past = 2;
    param.m = 5;
    param.linesearch = LBFGSpp::LBFGS_LINESEARCH_BACKTRACKING_WOLFE;
    
    // Create solver and function object
    LBFGSpp::LBFGSSolver<DataType> m_solver(param);
    
    double fx = 0.0;
    Eigen::VectorXd qNew = P*qDot;
    qNew.setZero();
   
    m_solver.minimize(objective, qNew, fx);
    mapStateEigen<1>(world) = P.transpose()*qNew;
    mapStateEigen<0>(world) = q + dt*P.transpose()*qNew;
    
    
}

template<typename DataType, typename MatrixAssembler, typename VectorAssembler>
using TimeStepperEulerImplicitBFGS = TimeStepper<DataType, TimeStepperImplEulerImplicitBFGS<DataType, MatrixAssembler, VectorAssembler> >;


#endif /* TimeStepperEulerImplicit_h */
