//
//  TimeStepperEulerImplicit.h
//  Gauss
//
//  Created by David Levin on 11/23/17.
//
//

#ifndef TimeStepperEulerImplicit_h
#define TimeStepperEulerImplicit_h

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

//TODO Solver Interface
namespace Gauss {

    //Given Initial state, step forward in time using linearly implicit Euler Integrator
    template<typename DataType, typename MatrixAssembler, typename VectorAssembler>
    class TimeStepperImplEulerImplicit
    {
    public:

        TimeStepperImplEulerImplicit(unsigned int num_iterations = 10000)  {
            m_num_iterations = num_iterations;
        }

        TimeStepperImplEulerImplicit(const TimeStepperImplEulerImplicit &toCopy) {

        }

        ~TimeStepperImplEulerImplicit() { }

        //Methods
        //init() //initial conditions will be set at the begining
        template<typename World>
        void step(World &world, double dt, double t);

        inline typename VectorAssembler::MatrixType & getLagrangeMultipliers() { return m_lagrangeMultipliers; }

    protected:

        MatrixAssembler m_massMatrix;
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
void TimeStepperImplEulerImplicit<DataType, MatrixAssembler, VectorAssembler>::step(World &world, double dt, double t) {

    //First two lines work around the fact that C++11 lambda can't directly capture a member variable.
    MatrixAssembler &massMatrix = m_massMatrix;
    MatrixAssembler &stiffnessMatrix = m_stiffnessMatrix;
    MatrixAssembler &AeqMatrix = m_Aeq;
    
    VectorAssembler &bVector = m_b;
    VectorAssembler &forceVector = m_forceVector;

    //Grab the state
    Eigen::VectorXd q = mapStateEigen<0>(world);
    Eigen::VectorXd qDot = mapStateEigen<1>(world);


    ///get mass matrix
    ASSEMBLEMATINIT(massMatrix, world.getNumQDotDOFs(), world.getNumQDotDOFs());
    ASSEMBLELIST(massMatrix, world.getSystemList(), getMassMatrix);
    ASSEMBLEEND(massMatrix);

    //we're going to build equality constraints into our gradient and hessian calcuations
    auto E = [&world, &massMatrix, &qDot](auto &a) { return (getEnergy(world) -
                                                             mapStateEigen<1>(world).transpose()*(*massMatrix)*qDot); };
    
    auto H = [&world, &massMatrix, &stiffnessMatrix, &dt, &qDot](auto &a)->auto & {
        //get stiffness matrix
        ASSEMBLEMATINIT(stiffnessMatrix, world.getNumQDotDOFs(), world.getNumQDotDOFs());
        ASSEMBLELIST(stiffnessMatrix, world.getSystemList(), getStiffnessMatrix);
        ASSEMBLELIST(stiffnessMatrix, world.getForceList(), getStiffnessMatrix);
        ASSEMBLEEND(stiffnessMatrix);
        
        (*stiffnessMatrix) *= -(dt*dt);
        (*stiffnessMatrix) += (*massMatrix);
        return stiffnessMatrix;
    };
    
    auto Aeq = [&world, &massMatrix, &AeqMatrix, &dt, &qDot](auto &a)->auto & {
        
        ASSEMBLEMATINIT(AeqMatrix, world.getNumConstraints(), world.getNumQDotDOFs());
        ASSEMBLELISTCONSTRAINT(AeqMatrix, world.getConstraintList(), getGradient);
        ASSEMBLEEND(AeqMatrix);
        
        return AeqMatrix;
    };

    auto g = [&world, &massMatrix, &forceVector, &dt, &qDot](auto &a) -> auto & {
        ASSEMBLEVECINIT(forceVector, world.getNumQDotDOFs());
        ASSEMBLELIST(forceVector, world.getForceList(), getForce);
        ASSEMBLELIST(forceVector, world.getSystemList(), getForce);
        ASSEMBLEEND(forceVector);
        
        (*forceVector).head(world.getNumQDotDOFs()) *= -dt;
        (*forceVector).head(world.getNumQDotDOFs()) += (*massMatrix)*(mapStateEigen<1>(world)-qDot);
        
        return forceVector;
    };
    
    auto b = [&world, &massMatrix, &bVector, &dt, &qDot](auto &a) -> auto & {
        ASSEMBLEVECINIT(bVector, world.getNumConstraints());
        ASSEMBLELISTCONSTRAINT(bVector, world.getConstraintList(), getDbDt);
        ASSEMBLEEND(bVector);
        
        return bVector;
    };
    
    auto update = [&world, &q, &qDot, &dt, &massMatrix](auto &dx) {

        //std::cout<<"NORM: "<<dx.head(world.getNumQDOFs()).norm()<<"\n";
        mapStateEigen<1>(world) = dx.head(world.getNumQDOFs());
        mapStateEigen<0>(world) = q + dt*mapStateEigen<1>(world);

    };

    //solve this using newton's method
    auto x0 = Eigen::VectorXd(world.getNumQDotDOFs()+world.getNumConstraints());
    x0.setZero();
    m_newton(x0, E, g, H, b, Aeq, update, 1e-4, m_num_iterations);

}

template<typename DataType, typename MatrixAssembler, typename VectorAssembler>
using TimeStepperEulerImplicit = TimeStepper<DataType, TimeStepperImplEulerImplicit<DataType, MatrixAssembler, VectorAssembler> >;

// #endif
#endif /* TimeStepperEulerImplicit_h */
