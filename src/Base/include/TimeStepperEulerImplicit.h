//
//  TimeStepperEulerImplicit.h
//  Gauss
//
//  Created by David Levin on 11/23/17.
//
//

#ifndef TimeStepperEulerImplicit_h
#define TimeStepperEulerImplicit_h

#ifdef GAUSS_PARDISO

#include <World.h>
#include <Assembler.h>
#include <TimeStepper.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <UtilitiesEigen.h>
#include <UtilitiesMATLAB.h>
#include <Eigen/SparseCholesky>
#include <SolverPardiso.h>
#include <Newton.h>
//TODO Solver Interface
namespace Gauss {
    
    //Given Initial state, step forward in time using linearly implicit Euler Integrator
    template<typename DataType, typename MatrixAssembler, typename VectorAssembler>
    class TimeStepperImplEulerImplicit
    {
    public:
        
        TimeStepperImplEulerImplicit()  {
            
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
        VectorAssembler m_forceVector;
        
        typename VectorAssembler::MatrixType m_lagrangeMultipliers;
        
        SolverPardiso<Eigen::SparseMatrix<DataType, Eigen::RowMajor> > m_pardiso;
        
    private:
    };
}

template<typename DataType, typename MatrixAssembler, typename VectorAssembler>
template<typename World>
void TimeStepperImplEulerImplicit<DataType, MatrixAssembler, VectorAssembler>::step(World &world, double dt, double t) {
    
    //First two lines work around the fact that C++11 lambda can't directly capture a member variable.
    MatrixAssembler &massMatrix = m_massMatrix;
    MatrixAssembler &stiffnessMatrix = m_stiffnessMatrix;
    VectorAssembler &forceVector = m_forceVector;
    SolverPardiso<Eigen::SparseMatrix<DataType, Eigen::RowMajor> > &solver = m_pardiso;
    
    //Grab the state
    Eigen::VectorXd q = mapStateEigen<0>(world);
    Eigen::VectorXd qDot = mapStateEigen<1>(world);
    
    
    ///get mass matrix
    ASSEMBLEMATINIT(massMatrix, world.getNumQDotDOFs()+world.getNumConstraints(), world.getNumQDotDOFs()+world.getNumConstraints());
    ASSEMBLELIST(massMatrix, world.getSystemList(), getMassMatrix);
    
    //add in the constraints
    ASSEMBLELISTOFFSET(massMatrix, world.getConstraintList(), getGradient, world.getNumQDotDOFs(), 0);
    ASSEMBLELISTOFFSETTRANSPOSE(massMatrix, world.getConstraintList(), getGradient, 0, world.getNumQDotDOFs());
    ASSEMBLEEND(massMatrix);
    
    //std::cout<<"F: "<<(*forceVector)<<"\n";
    //setup RHS

    
    //toMatlab(*massMatrix, "./testMHex.txt");
    //toMatlab(*stiffnessMatrix, "./testKHex.txt");
    //toMatlab(*forceVector, "./testfHex.txt");
    
    //std::cout<<"F: \n"<<(*forceVector)<<"\n";
    
    
    //Eigen::SparseMatrix<DataType, Eigen::RowMajor> systemMatrix = (*m_massMatrix)- dt*dt*(*m_stiffnessMatrix);
    
    //we're going to build equality constraints into our gradient and hessian calcuations
    auto E = [&world](auto &a) { return getEnergy(world); };
    
    auto H = [&world, &massMatrix, &stiffnessMatrix, &dt, &qDot](auto &a)->auto & {
        //get stiffness matrix
        ASSEMBLEMATINIT(stiffnessMatrix, world.getNumQDotDOFs()+world.getNumConstraints(), world.getNumQDotDOFs()+world.getNumConstraints());
        ASSEMBLELIST(stiffnessMatrix, world.getSystemList(), getStiffnessMatrix);
        ASSEMBLELIST(stiffnessMatrix, world.getForceList(), getStiffnessMatrix);
        
        ASSEMBLEEND(stiffnessMatrix);
        
        (*stiffnessMatrix) *= -(dt*dt);
        (*stiffnessMatrix) += (*massMatrix);
        return (*stiffnessMatrix);
    };
    
    auto g = [&world, &massMatrix, &forceVector, &dt, &qDot](auto &a) -> auto & {
        ASSEMBLEVECINIT(forceVector, world.getNumQDotDOFs()+world.getNumConstraints());
        ASSEMBLELIST(forceVector, world.getForceList(), getForce);
        ASSEMBLELIST(forceVector, world.getSystemList(), getForce);
        ASSEMBLEEND(forceVector);
        
        (*forceVector).head(world.getNumQDotDOFs()) *= -dt;
        (*forceVector).head(world.getNumQDotDOFs()) += (*massMatrix).block(0,0,world.getNumQDotDOFs(),world.getNumQDotDOFs())*(mapStateEigen<1>(world)-qDot);
        
        return (*forceVector);
    };
    
    auto solve = [&solver](auto &A, auto &b)->auto & {
        return solver.solve(A,b);
    };
    
    auto update = [&world, &q, &qDot, &dt, &massMatrix](auto &dx) {
        
        //std::cout<<"NORM: "<<dx.head(world.getNumQDOFs()).norm()<<"\n";
        mapStateEigen<1>(world) = dx.head(world.getNumQDOFs());
        mapStateEigen<0>(world) = q + dt*mapStateEigen<1>(world);
    
    };
    
    //solve this using newton's method
    auto x0 = Eigen::VectorXd(world.getNumQDotDOFs()+world.getNumConstraints());
    x0.setZero();
    mapStateEigen<1>(world).setZero();
    //x0.head(world.getNumQDotDOFs()) = mapStateEigen<1>(world);
    Optimization::newton(E, g, H, solve, x0, update, 1e-4, 10000);
    //std::cout<<mapStateEigen<1>(world)<<"\n";
    
    
}

template<typename DataType, typename MatrixAssembler, typename VectorAssembler>
using TimeStepperEulerImplicit = TimeStepper<DataType, TimeStepperImplEulerImplicit<DataType, MatrixAssembler, VectorAssembler> >;

#endif
#endif /* TimeStepperEulerImplicit_h */
