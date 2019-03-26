//
//  TimeStepperEulerImplicit.h
//  Gauss
//
//  Created by David Levin on 11/23/17.
//
//

#ifndef TimeStepperStatic_h
#define TimeStepperStatic_h

// #ifdef GAUSS_PARDISO

#include <World.h>
#include <Assembler.h>
#include <TimeStepper.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <UtilitiesEigen.h>
#include <UtilitiesMATLAB.h>
#include <Eigen/SparseCholesky>
#include <SolverPardiso.h>
#include <GaussOptimizationAdapters.h>

//This takes a single static time step by just minimizing the potential energy of the object
namespace Gauss {
    
    //Given Initial state, step forward in time using linearly implicit Euler Integrator
    template<typename DataType, typename MatrixAssembler, typename VectorAssembler>
    class TimeStepperImplStatic
    {
    public:
        
        TimeStepperImplStatic(unsigned int num_iterations = 1000)  {
            m_num_iterations = num_iterations;
        }
        
        TimeStepperImplStatic(const TimeStepperImplStatic &toCopy) {
            
        }
        
        ~TimeStepperImplStatic() { }
        
        //Methods
        //init() //initial conditions will be set at the begining
        template<typename World>
        void step(World &world, double dt, double t);
        
        inline typename VectorAssembler::MatrixType & getLagrangeMultipliers() { return m_lagrangeMultipliers; }
        
    protected:
        
        MatrixAssembler m_stiffnessMatrix;
        MatrixAssembler m_Aeq;
        VectorAssembler m_forceVector;
        VectorAssembler m_b;
        
        unsigned int m_num_iterations;
        typename VectorAssembler::MatrixType m_lagrangeMultipliers;
        Optimization::NewtonSearchWithBackTracking<DataType> m_newton;
        
    private:
    };
}

template<typename DataType, typename MatrixAssembler, typename VectorAssembler>
template<typename World>
void TimeStepperImplStatic<DataType, MatrixAssembler, VectorAssembler>::step(World &world, double dt, double t) {
    
    //First two lines work around the fact that C++11 lambda can't directly capture a member variable.
    MatrixAssembler &stiffnessMatrix = m_stiffnessMatrix;
    MatrixAssembler &AeqMatrix = m_Aeq;
    
    VectorAssembler &bVector = m_b;
    VectorAssembler &forceVector = m_forceVector;
    
    //Grab the state
    Eigen::VectorXd q = mapStateEigen<0>(world);
    
    //we're going to build equality constraints into our gradient and hessian calcuations
    auto E = [&world](auto &a) { //return (getEnergy(world) -
        return getEnergy(world);
    };
    
    auto H = [&world, &stiffnessMatrix](auto &a)->auto & {
        //get stiffness matrix
        ASSEMBLEMATINIT(stiffnessMatrix, world.getNumQDotDOFs(), world.getNumQDotDOFs());
        ASSEMBLELIST(stiffnessMatrix, world.getSystemList(), getStiffnessMatrix);
        ASSEMBLELIST(stiffnessMatrix, world.getForceList(), getStiffnessMatrix);
        ASSEMBLEEND(stiffnessMatrix);
        
        (*stiffnessMatrix) *= -1.0;
        
        return stiffnessMatrix;
    };
    
    auto Aeq = [&world, &AeqMatrix](auto &a)->auto & {
        
        ASSEMBLEMATINIT(AeqMatrix, world.getNumConstraints(), world.getNumQDotDOFs());
        ASSEMBLELISTCONSTRAINT(AeqMatrix, world.getConstraintList(), getGradient);
        ASSEMBLEEND(AeqMatrix);
        
        return AeqMatrix;
    };
    
    auto g = [&world, &forceVector](auto &a) -> auto & {
        ASSEMBLEVECINIT(forceVector, world.getNumQDotDOFs());
        ASSEMBLELIST(forceVector, world.getForceList(), getForce);
        ASSEMBLELIST(forceVector, world.getSystemList(), getForce);
        ASSEMBLEEND(forceVector);
        
        (*forceVector).head(world.getNumQDotDOFs()) *= -1.0;
        
        return forceVector;
    };
    
    auto b = [&world, &bVector](auto &a) -> auto & {
        ASSEMBLEVECINIT(bVector, world.getNumConstraints());
        ASSEMBLELISTCONSTRAINT(bVector, world.getConstraintList(), getDbDt);
        ASSEMBLEEND(bVector);
        
        return bVector;
    };
    
    auto update = [&world, &q](auto &dx) {
        
        mapStateEigen<0>(world) = dx.head(world.getNumQDOFs());
        
    };
    
    //solve this using newton's method
    auto x0 = Eigen::VectorXd(world.getNumQDOFs()+world.getNumConstraints());
    x0.setZero();
    m_newton(x0, E, g, H, b, Aeq, update, 1e-4, m_num_iterations);
    
    
}

template<typename DataType, typename MatrixAssembler, typename VectorAssembler>
using TimeStepperStatic = TimeStepper<DataType, TimeStepperImplStatic<DataType, MatrixAssembler, VectorAssembler> >;

// #endif
#endif /* TimeStepperStatic_h */

