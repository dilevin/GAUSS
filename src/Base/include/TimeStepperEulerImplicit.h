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

        TimeStepperImplEulerImplicit(unsigned int num_iterations)  {
            m_num_iterations = num_iterations;
#ifdef GAUSS_PARDISO
            std::cout<<"TimeStepperEuler.h will use Pardiso solver"<<std::endl;
#else
            std::cout<<"TimeStepperEuler.h will use Eigen's solver"<<std::endl;
#endif
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
#ifdef GAUSS_PARDISO
        SolverPardiso<Eigen::SparseMatrix<DataType, Eigen::RowMajor> > m_pardiso;
#else
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > m_eigensolver;
#endif

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
#ifdef GAUSS_PARDISO
    SolverPardiso<Eigen::SparseMatrix<DataType, Eigen::RowMajor> > &solver = m_pardiso;
#else
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > &solver = m_eigensolver;
    Eigen::VectorXd eigenSolverResult;
#endif
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
    auto E = [&world, &massMatrix, &qDot](auto &a) { return (getEnergy(world) -
                                                             mapStateEigen<1>(world).transpose()*(*massMatrix).block(0,0,world.getNumQDotDOFs(),world.getNumQDotDOFs())*qDot); };

    /*auto H = [&world, &massMatrix, &stiffnessMatrix, &dt, &qDot](auto &a)->auto & {
        //get stiffness matrix
        ASSEMBLEMATINIT(stiffnessMatrix, world.getNumQDotDOFs()+world.getNumConstraints(), world.getNumQDotDOFs()+world.getNumConstraints());
        ASSEMBLELIST(stiffnessMatrix, world.getSystemList(), getStiffnessMatrix);
        ASSEMBLELIST(stiffnessMatrix, world.getForceList(), getStiffnessMatrix);

        ASSEMBLEEND(stiffnessMatrix);

        (*stiffnessMatrix) *= -(dt*dt);
        (*stiffnessMatrix) += (*massMatrix);
        return (*stiffnessMatrix);
    };*/
    
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
        //get stiffness matrix
        ASSEMBLEMATINIT(AeqMatrix, world.getNumConstraints(), world.getNumQDotDOFs());
        ASSEMBLELIST(AeqMatrix, world.getConstraintList(), getGradient);
        ASSEMBLEEND(AeqMatrix);
        
        return AeqMatrix;
    };

    /*auto g = [&world, &massMatrix, &forceVector, &dt, &qDot](auto &a) -> auto & {
        ASSEMBLEVECINIT(forceVector, world.getNumQDotDOFs()+world.getNumConstraints());
        ASSEMBLELIST(forceVector, world.getForceList(), getForce);
        ASSEMBLELIST(forceVector, world.getSystemList(), getForce);
        ASSEMBLELISTOFFSET(forceVector, world.getConstraintList(), getDbDt, world.getNumQDotDOFs(), 0);
        ASSEMBLEEND(forceVector);

        (*forceVector).head(world.getNumQDotDOFs()) *= -dt;
        (*forceVector).tail(world.getNumConstraints()) *= -1.0;
        (*forceVector).head(world.getNumQDotDOFs()) += (*massMatrix).block(0,0,world.getNumQDotDOFs(),world.getNumQDotDOFs())*(mapStateEigen<1>(world)-qDot);

        return (*forceVector);
    };*/
    
    auto g = [&world, &massMatrix, &forceVector, &dt, &qDot](auto &a) -> auto & {
        ASSEMBLEVECINIT(forceVector, world.getNumQDotDOFs());
        ASSEMBLELIST(forceVector, world.getForceList(), getForce);
        ASSEMBLELIST(forceVector, world.getSystemList(), getForce);
        ASSEMBLEEND(forceVector);
        
        (*forceVector).head(world.getNumQDotDOFs()) *= -dt;
        (*forceVector).head(world.getNumQDotDOFs()) += (*massMatrix).block(0,0,world.getNumQDotDOFs(),world.getNumQDotDOFs())*(mapStateEigen<1>(world)-qDot);
        
        return forceVector;
    };
    
    auto b = [&world, &massMatrix, &bVector, &dt, &qDot](auto &a) -> auto & {
        ASSEMBLEVECINIT(bVector, world.getNumConstraints());
        ASSEMBLELIST(bVector, world.getConstraintList, getDbDt);
        ASSEMBLEEND(bVector);
        
        return bVector;
    };
    
#ifdef GAUSS_PARDISO

    auto solve = [&solver](auto &A, auto &b)->auto & {
        return solver.solve(A,b);
    };
#else
    auto solve = [&solver, &eigenSolverResult](auto &A, auto &b)->auto & {
        solver.compute(A);
        if(solver.info()!=Eigen::Success) {
            // decomposition failed
            // assert(1 == 0);
            std::cout<<"Decomposition Failed \n";
            exit(1);
        }

        if(solver.info()!=Eigen::Success) {
            // solving failed
            // assert(1 == 0);
            std::cout<<"Solve Failed \n";
            exit(1);
        }
        eigenSolverResult = solver.solve(b);
        return eigenSolverResult;
    };
#endif

    auto update = [&world, &q, &qDot, &dt, &massMatrix](auto &dx) {

        //std::cout<<"NORM: "<<dx.head(world.getNumQDOFs()).norm()<<"\n";
        mapStateEigen<1>(world) = dx.head(world.getNumQDOFs());
        mapStateEigen<0>(world) = q + dt*mapStateEigen<1>(world);

    };

    //solve this using newton's method
    auto x0 = Eigen::VectorXd(world.getNumQDotDOFs()+world.getNumConstraints());
    x0.setZero();
    //x0.head(world.getNumQDotDOFs()) = mapStateEigen<1>(world);
    //Optimization::newton(E, g, H, solve, x0, update, 1e-4, m_num_iterations);
    //std::cout<<mapStateEigen<1>(world)<<"\n";
    Optimization::minimizeWorldNewton(x0, E, g, H, b, Aeq, update, 1e-4, m_num_iterations);


}

template<typename DataType, typename MatrixAssembler, typename VectorAssembler>
using TimeStepperEulerImplicit = TimeStepper<DataType, TimeStepperImplEulerImplicit<DataType, MatrixAssembler, VectorAssembler> >;

// #endif
#endif /* TimeStepperEulerImplicit_h */
