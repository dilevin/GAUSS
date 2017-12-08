//
//  TimeStepperEulerImplicitLinearCollisions.h
//  Gauss
//
//  Created by David Levin on 12/6/17.
//
//

#ifndef TimeStepperEulerImplicitLinearCollisions_h
#define TimeStepperEulerImplicitLinearCollisions_h

#include <CollisionDetector.h>
#include <CollisionsFloor.h>
#include <World.h>
#include <Assembler.h>
#include <TimeStepper.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <UtilitiesEigen.h>
#include <UtilitiesMATLAB.h>
#include <Eigen/SparseCholesky>
#include <SolverPardiso.h>

//solver
#include <igl/active_set.h>

namespace Gauss {
    namespace Collisions {
        //Given Initial state, step forward in time using linearly implicit Euler Integrator
        template<typename DataType, typename MatrixAssembler, typename VectorAssembler>
        class TimeStepperImplEulerImplicitLinearCollisions
        {
        public:
            
            TimeStepperImplEulerImplicitLinearCollisions() {
                
            }
            
            TimeStepperImplEulerImplicitLinearCollisions(const TimeStepperImplEulerImplicitLinearCollisions &toCopy) {
                
            }
            
            ~TimeStepperImplEulerImplicitLinearCollisions() { }
            
            //Methods
            //init() //initial conditions will be set at the begining
            template<typename World>
            void step(World &world, double dt, double t);
            
            inline typename VectorAssembler::MatrixType & getLagrangeMultipliers() { return m_lagrangeMultipliers; }
            
        protected:
            
            MatrixAssembler m_massMatrix;
            MatrixAssembler m_stiffnessMatrix;
            MatrixAssembler m_collisionConstraints;
            VectorAssembler m_forceVector;
            
            //storage for lagrange multipliers
            typename VectorAssembler::MatrixType m_lagrangeMultipliers;
            
#ifdef GAUSS_PARDISO
            
            SolverPardiso<Eigen::SparseMatrix<DataType, Eigen::RowMajor> > m_pardiso;
            
#endif
            
        private:
        };
    
    template<typename DataType, typename MatrixAssembler, typename VectorAssembler>
    template<typename World>
    void TimeStepperImplEulerImplicitLinearCollisions<DataType, MatrixAssembler, VectorAssembler>::step(World &world, double dt, double t) {
        
        //First two lines work around the fact that C++11 lambda can't directly capture a member variable.
        MatrixAssembler &massMatrix = m_massMatrix;
        MatrixAssembler &collisions = m_collisionConstraints;
        MatrixAssembler &stiffnessMatrix = m_stiffnessMatrix;
        VectorAssembler &forceVector = m_forceVector;
        
        //get mass matrix
        ASSEMBLEMATINIT(massMatrix, world.getNumQDotDOFs(), world.getNumQDotDOFs()+world.getNumConstraints());
        ASSEMBLELIST(massMatrix, world.getSystemList(), getMassMatrix);
        ASSEMBLEEND(massMatrix);
        
        
        //get stiffness matrix
        ASSEMBLEMATINIT(stiffnessMatrix, world.getNumQDotDOFs(), world.getNumQDotDOFs());
        ASSEMBLELIST(stiffnessMatrix, world.getSystemList(), getStiffnessMatrix);
        ASSEMBLELIST(stiffnessMatrix, world.getForceList(), getStiffnessMatrix);
        ASSEMBLEEND(stiffnessMatrix);
        
        ASSEMBLEVECINIT(forceVector, world.getNumQDotDOFs());
        ASSEMBLELIST(forceVector, world.getForceList(), getForce);
        ASSEMBLELIST(forceVector, world.getSystemList(), getForce);
        ASSEMBLEEND(forceVector);
        
        //collision constraints
        
        //make sure collisions have been detected by running the constraint update function
        world.updateInequalityConstraints();
        ASSEMBLEMATINIT(collisions, world.getNumInequalityConstraints(), world.getNumQDotDOFs());
        ASSEMBLELISTCONSTRAINT(collisions, world.getInequalityConstraintList(), getGradient);
        ASSEMBLEEND(collisions);
    
        
        //Grab the state
        Eigen::Map<Eigen::VectorXd> q = mapStateEigen<0>(world);
        Eigen::Map<Eigen::VectorXd> qDot = mapStateEigen<1>(world);
        
        //setup RHS
        (*forceVector) = -((*massMatrix)*qDot + dt*(*forceVector));
        
        Eigen::VectorXd x0;
        Eigen::SparseMatrix<DataType> systemMatrix = (*m_massMatrix)- dt*dt*(*m_stiffnessMatrix);
        
        
        //solve using libigl active set solver
        Eigen::VectorXi known(0);
        Eigen::VectorXd bKnown(0);
        Eigen::SparseMatrix<DataType> Aeq;
        Eigen::VectorXd beq;
        Eigen::VectorXd b;
        Eigen::SparseMatrix<DataType> Aineq = -(*collisions);
        b.resize((*collisions).rows(),1);
        b.setZero();
        Eigen::VectorXd lx;
        Eigen::VectorXd ux;
        igl::active_set_params params;
        
        Eigen::VectorXd qDotTmp = qDot;
        active_set(systemMatrix, (*forceVector), known, bKnown, Aeq, beq, Aineq, b, lx, ux, params, qDotTmp);
        qDot = qDotTmp;
        q = q + dt*qDot;
        
        /*active_set(
                   const Eigen::SparseMatrix<AT>& A,
                   const Eigen::PlainObjectBase<DerivedB> & B,
                   const Eigen::PlainObjectBase<Derivedknown> & known,
                   const Eigen::PlainObjectBase<DerivedY> & Y,
                   const Eigen::SparseMatrix<AeqT>& Aeq,
                   const Eigen::PlainObjectBase<DerivedBeq> & Beq,
                   const Eigen::SparseMatrix<AieqT>& Aieq,
                   const Eigen::PlainObjectBase<DerivedBieq> & Bieq,
                   const Eigen::PlainObjectBase<Derivedlx> & lx,
                   const Eigen::PlainObjectBase<Derivedux> & ux,
                   const igl::active_set_params & params,
                   Eigen::PlainObjectBase<DerivedZ> & Z
                   );*.
        
         
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
    }*/
    }
    
    template<typename DataType, typename MatrixAssembler, typename VectorAssembler>
    using TimeStepperEulerImplictLinearCollisions = TimeStepper<DataType, TimeStepperImplEulerImplicitLinearCollisions<DataType, MatrixAssembler, VectorAssembler> >;
    }
}
#endif /* TimeStepperEulerImplicitLinearCollisions_h */
