#ifndef TimeStepperNewmarkLinear_h
#define TimeStepperNewmarkLinear_h

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


namespace Gauss {
    
    template<typename DataType, typename MatrixAssembler, typename VectorAssembler>
    class TimeStepperImpNewmark {
    public:
        
        TimeStepperImpNewmark(unsigned int num_iterations = 100) {
            //m_P = P;
            m_num_iterations = num_iterations;
            
        }
        
        TimeStepperImpNewmark(const TimeStepperImpNewmark &toCopy) {
            
        }
        
        ~TimeStepperImpNewmark() {
        }
        
        //Methods
        template<typename World>
        void step(World &world, double dt, double t);
        
        inline typename VectorAssembler::MatrixType & getLagrangeMultipliers() { return m_lagrangeMultipliers; }
        
    protected:
        
        bool initialized = false;
        
        unsigned int m_num_iterations;
        MatrixAssembler m_massMatrix;
        MatrixAssembler m_stiffnessMatrix;
        MatrixAssembler m_Aeq;
        VectorAssembler m_forceVector;
        VectorAssembler m_b;
        
        //Eigen::SparseMatrix<double> m_P;
        
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
        
        //Newton solver 
        Optimization::NewtonSearchWithBackTracking<DataType> m_newton;
        //Optimization::NewtonSearchWithAndersonAcceleration<DataType> m_newton;
        
        //storage for lagrange multipliers
        typename VectorAssembler::MatrixType m_lagrangeMultipliers;
        
    };
}

template<typename DataType, typename MatrixAssembler, typename VectorAssembler>
template<typename World>
void TimeStepperImpNewmark<DataType, MatrixAssembler, VectorAssembler>::step(World &world, double dt, double t) {
    
    //First two lines work around the fact that C++11 lambda can't directly capture a member variable.
    Eigen::VectorXx<DataType> delta;
    MatrixAssembler &massMatrix = m_massMatrix;
    MatrixAssembler &stiffnessMatrix = m_stiffnessMatrix;
    MatrixAssembler &AeqMatrix = m_Aeq;
    
    VectorAssembler &bVector = m_b;
    VectorAssembler &forceVector = m_forceVector;
    
    // precompute and prefactor Hessian
    if (!initialized) {
        //assume mass matrix is constant
        getMassMatrix(m_massMatrix, world);
        initialized = true;
    }
    
    // assemble b vector
    Eigen::VectorXd q = mapStateEigen<0>(world);
    Eigen::VectorXd qDot = mapStateEigen<1>(world);
    AssemblerEigenVector<double> force;
    getForceVector(force, world);
    Eigen::VectorXx<DataType>  b = dt * (*m_massMatrix)*qDot + (dt*dt / 4.0)*(*force);
    
    //initialize delta
    delta.resize(world.getNumQDotDOFs()+world.getNumConstraints(),1);
    delta.setZero();
    
    //lambdas for newton solver
    //we're going to build equality constraints into our gradient and hessian calcuations
    auto E = [&world, &massMatrix, &b, &delta, &dt](auto &a) { return (0.25*dt*dt)*(getStrainEnergy(world) + getBodyForceEnergy(world)) - a.head(world.getNumQDOFs()).dot(b) + 0.5*a.head(world.getNumQDOFs()).transpose()*(*massMatrix)*a.head(world.getNumQDOFs()); };
    
    auto H = [&world, &massMatrix, &stiffnessMatrix, &dt, &qDot](auto &a)->auto & {
        //get stiffness matrix
        ASSEMBLEMATINIT(stiffnessMatrix, world.getNumQDotDOFs(), world.getNumQDotDOFs());
        ASSEMBLELIST(stiffnessMatrix, world.getSystemList(), getStiffnessMatrix);
        ASSEMBLELIST(stiffnessMatrix, world.getForceList(), getStiffnessMatrix);
        ASSEMBLEEND(stiffnessMatrix);
        
        (*stiffnessMatrix) *= -0.25*(dt*dt);
        (*stiffnessMatrix) += (*massMatrix);
        return stiffnessMatrix;
    };
    
    auto Aeq = [&world, &massMatrix, &AeqMatrix, &dt, &qDot](auto &a)->auto & {
        ASSEMBLEMATINIT(AeqMatrix, world.getNumConstraints(), world.getNumQDotDOFs());
        ASSEMBLELISTCONSTRAINT(AeqMatrix, world.getConstraintList(), getGradient);
        ASSEMBLEEND(AeqMatrix);
        
        return AeqMatrix;
    };
    
    auto g = [&world, &massMatrix, &forceVector, &b, &delta, &dt](auto &a) -> auto & {
        //get force vector at current configuration (qt + delta)
        //configuration is set in the update method below via the callback function sent to the newton's solver
        ASSEMBLEVECINIT(forceVector, world.getNumQDotDOFs());
        ASSEMBLELIST(forceVector, world.getForceList(), getForce);
        ASSEMBLELIST(forceVector, world.getSystemList(), getForce);
        ASSEMBLEEND(forceVector);
        
        (*forceVector).head(world.getNumQDotDOFs()) *= -0.25*(dt*dt);
        (*forceVector).head(world.getNumQDotDOFs()) += ((*massMatrix)*a.head(world.getNumQDOFs()) - b);
        
        return forceVector;
    };
    
    auto beq = [&world, &massMatrix, &bVector, &dt, &qDot](auto &a) -> auto & {
        ASSEMBLEVECINIT(bVector, world.getNumConstraints());
        ASSEMBLELISTCONSTRAINT(bVector, world.getConstraintList(), getDbDt);
        ASSEMBLEEND(bVector);
        
        return bVector;
    };
    
    auto update = [&world, &q, &qDot, &delta, &dt](auto &dx) {
        
        mapStateEigen<0>(world) = q + dx.head(world.getNumQDOFs());
        mapStateEigen<1>(world) = (2.0/dt)*dx.head(world.getNumQDOFs()) - qDot;
    };
    
    //solve for delta
    m_newton(delta, E, g, H, beq, Aeq, update, 1e-4, m_num_iterations);
    
    // update state
    q = q + delta.head(world.getNumQDOFs());
    qDot = (2.0 / dt)*delta.head(world.getNumQDOFs()) - qDot;
}

template<typename DataType, typename MatrixAssembler, typename VectorAssembler>
using TimeStepperNewmark = TimeStepper<DataType, TimeStepperImpNewmark<DataType, MatrixAssembler, VectorAssembler> >;

#endif //TimeStepperNewmarkLinear_h

