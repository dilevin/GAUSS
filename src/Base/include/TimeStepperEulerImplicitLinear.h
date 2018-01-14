//
//  TimeStepperEulerImplicitLinear.h
//  Gauss
//
//  Created by David Levin on 2/10/17.
//
//

#ifndef TimeStepperEulerImplicitLinear_h
#define TimeStepperEulerImplicitLinear_h

#include <World.h>
#include <Assembler.h>
#include <TimeStepper.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <UtilitiesEigen.h>
#include <UtilitiesMATLAB.h>
#include <Eigen/SparseCholesky>
#include <SolverPardiso.h>

//TODO Solver Interface
namespace Gauss {
    
    //Given Initial state, step forward in time using linearly implicit Euler Integrator
    template<typename DataType, typename MatrixAssembler, typename VectorAssembler>
    class TimeStepperImplEulerImplicitLinear
    {
    public:
        
        TimeStepperImplEulerImplicitLinear() {
            
        }
        
        TimeStepperImplEulerImplicitLinear(const TimeStepperImplEulerImplicitLinear &toCopy) {
            
        }
        
        ~TimeStepperImplEulerImplicitLinear() { }
        
        //Methods
        //init() //initial conditions will be set at the begining
        template<typename World>
        void step(World &world, double dt, double t);
        
        inline typename VectorAssembler::MatrixType & getLagrangeMultipliers() { return m_lagrangeMultipliers; }
        
    protected:
        
        MatrixAssembler m_massMatrix;
        MatrixAssembler m_stiffnessMatrix;
        VectorAssembler m_forceVector;
        
        //storage for lagrange multipliers
        typename VectorAssembler::MatrixType m_lagrangeMultipliers;
        
#ifdef GAUSS_PARDISO
        
        SolverPardiso<Eigen::SparseMatrix<DataType, Eigen::RowMajor> > m_pardiso;
        
#endif
        
    private:
    };
}

template<typename DataType, typename MatrixAssembler, typename VectorAssembler>
template<typename World>
void TimeStepperImplEulerImplicitLinear<DataType, MatrixAssembler, VectorAssembler>::step(World &world, double dt, double t) {

    //First two lines work around the fact that C++11 lambda can't directly capture a member variable.
    MatrixAssembler &massMatrix = m_massMatrix;
    MatrixAssembler &stiffnessMatrix = m_stiffnessMatrix;
    VectorAssembler &forceVector = m_forceVector;
    
    //get mass matrix
    ASSEMBLEMATINIT(massMatrix, world.getNumQDotDOFs()+world.getNumConstraints(), world.getNumQDotDOFs()+world.getNumConstraints());
    ASSEMBLELIST(massMatrix, world.getSystemList(), getMassMatrix);
    
    //add in the constraints
    ASSEMBLELISTOFFSET(massMatrix, world.getConstraintList(), getGradient, world.getNumQDotDOFs(), 0);
    ASSEMBLELISTOFFSETTRANSPOSE(massMatrix, world.getConstraintList(), getGradient, 0, world.getNumQDotDOFs());
    ASSEMBLEEND(massMatrix);
    
    
    //get stiffness matrix
    ASSEMBLEMATINIT(stiffnessMatrix, world.getNumQDotDOFs()+world.getNumConstraints(), world.getNumQDotDOFs()+world.getNumConstraints());
    ASSEMBLELIST(stiffnessMatrix, world.getSystemList(), getStiffnessMatrix);
    ASSEMBLELIST(stiffnessMatrix, world.getForceList(), getStiffnessMatrix);
    ASSEMBLEEND(stiffnessMatrix);
    
    ASSEMBLEVECINIT(forceVector, world.getNumQDotDOFs()+world.getNumConstraints());
    ASSEMBLELIST(forceVector, world.getForceList(), getForce);
    ASSEMBLELIST(forceVector, world.getSystemList(), getForce);
    ASSEMBLELISTOFFSET(forceVector, world.getConstraintList(), getFunction, world.getNumQDotDOFs(), 0);
    ASSEMBLEEND(forceVector);

    //toMatlab(*constraintMatrix, "./testConstraintMat");
    //std::cout<<*constraintVector<<"\n";
    //ASSEMBLEMAT(world, massMatrix, getNumQDotDOFs, getNumQDotDOFs, getMassMatrix);
    //ASSEMBLEMAT(world, stiffnessMatrix, getNumQDotDOFs, getNumQDotDOFs, getStiffnessMatrix);
   // ASSEMBLEMAT(world, constraintMatrix, getNumConstraints, getNumQDotDOFs, getGradient);
    
    //ASSEMBLEVEC(world, forceVector, getNumQDotDOFs, getForce);
    //ASSEMBLEVEC(world, constraintVector, getNumConstraints, getFunction)
    
    //Grab the state
    Eigen::Map<Eigen::VectorXd> q = mapStateEigen<0>(world);
    Eigen::Map<Eigen::VectorXd> qDot = mapStateEigen<1>(world);
    

    //std::cout<<"F: "<<(*forceVector)<<"\n";
    //setup RHS
    (*forceVector).head(world.getNumQDotDOFs()) = (*massMatrix).block(0,0, world.getNumQDotDOFs(), world.getNumQDotDOFs())*qDot + dt*(*forceVector).head(world.getNumQDotDOFs());
    
    //toMatlab(*massMatrix, "./testMHex.txt");
    //toMatlab(*stiffnessMatrix, "./testKHex.txt");
    //toMatlab(*forceVector, "./testfHex.txt");
    
    //std::cout<<"F: \n"<<(*forceVector)<<"\n";
    
    Eigen::VectorXd x0;
    Eigen::SparseMatrix<DataType, Eigen::RowMajor> systemMatrix = (*m_massMatrix)- dt*dt*(*m_stiffnessMatrix);
    
#ifdef GAUSS_PARDISO
    m_pardiso.symbolicFactorization(systemMatrix);
    m_pardiso.numericalFactorization();
    m_pardiso.solve(*forceVector);
    x0 = m_pardiso.getX();
    m_pardiso.cleanup();
#else
    //solve system (Need interface for solvers but for now just use Eigen LLt)
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
    solver.compute(systemMatrix);
    
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
    
    qDot = x0.head(world.getNumQDotDOFs());
    
    m_lagrangeMultipliers = x0.tail(world.getNumConstraints());
    
    //update state
    q = q + dt*qDot;
    
    //std::cout<<"Q: "<<q<<"\n";
}

template<typename DataType, typename MatrixAssembler, typename VectorAssembler>
using TimeStepperEulerImplicitLinear = TimeStepper<DataType, TimeStepperImplEulerImplicitLinear<DataType, MatrixAssembler, VectorAssembler> >;





#endif /* TimeStepperEulerImplicitLinear_h */
