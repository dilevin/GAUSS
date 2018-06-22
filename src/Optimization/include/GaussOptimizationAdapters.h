//
//  GaussOptimizationAdapters.h
//  Gauss
//
//  Created by David Levin on 6/13/18.
//

#ifndef GaussOptimizationAdapters_h
#define GaussOptimizationAdapters_h

#include <igl/cat.h>
#include <climits>
#include <Newton.h>
#include <UtilitiesEigen.h>

namespace Gauss {
    namespace Optimization {
        //functions that produce functors for plugging Gauss assemblers directly into my optimizers

        //Energy functor
        template<typename DataType, template<typename A> class Assembler, typename Func>
        class AssembledFunction {
            
        public:
            
            inline AssembledFunction(Func &func, Assembler<DataType> &assembler, unsigned int m, unsigned int n=1) : m_func(func), m_assembler(assembler), m_m(m), m_n(n) { }
            
            //here the input is the velocity we are solving for
            //copying matlabs interface where the function operator returns the energy value, the gradient and (in my code, the hessian approximation)
            //need more flexibility in my assemblings (scalar multiples etc)
            template<typename Vector>
            inline auto operator()(Vector &qDot, Assembler<DataType> &assembler) { //this takes the lambda as a parameter
                
                m_func(assembler);
                return assembler;
                
            }
        
            inline unsigned int getM() { return m_m; }
            
            inline unsigned int getN() {return m_n; }
            
            template<typename Vector>
            inline auto operator()(Vector &qDot) {
                
                Assembler<DataType> &assm = m_assembler;
                
                //init assembler
                ASSEMBLEMATINIT(m_assembler, m_m, m_n);
                this->operator()(qDot, assm);
                ASSEMBLEEND(m_assembler);
                
                return m_assembler;
                
            }
            
        protected:
            Assembler<DataType> &m_assembler;
            Func &m_func;
            unsigned int m_m, m_n; //assembler sizes
            
        };


        //Handy update rules for optimization
        //Velocity level dynamic update
        /*template<typename World, typename DataType>
        class UpdateQDot {
        public:
            UpdateQDot(World &world, DataType dt) : m_world(world), m_dt(dt) { }
            
            template<typename Vector>
            void operator()(Vector &qDot) {
                mapStateEigen<1>(m_world) = qDot.head(world.getNumQDOFs());
                mapStateEigen<0>(m_world) = q + qDot*mapStateEigen<1>(world);
            }
            
        protected:
            World &m_world;
            DataType m_dt;
            
        private:
        };*/
        
        //Newtons Step Direction from Assembler input (allows assembling everything ing place
        template<typename DataType>
        class DirectionNewtonAssembler {
          
        public:
            inline DirectionNewtonAssembler()
            #ifdef GAUSS_PARDISO
            : m_solver(11)
            #endif
            {
                
            }
            
            //takes in assembled matrices and returns the newton step direction
            template<typename Hessian, typename Gradient, typename Ceq, typename JacobianEq, typename Vector>
            inline decltype(auto) operator()(Hessian &H, Gradient &g, Ceq &ceq, JacobianEq &Jeq, Vector &x0) {
                
                
                
                //Use symmetric matrix format to make life easier
                unsigned int nnzH, nnzJ;
                nnzH = (*H).nonZeros();
                nnzJ = (*Jeq).outerIndexPtr()[(*Jeq).rows()];
                
                m_KKT.resize((*H).rows(), (*H).cols());
                m_KKT.reserve(nnzH + 2*nnzJ);
                m_KKT = (*H);
                
                //Build KKT using Eigen operations
                m_KKT.conservativeResize((*H).rows()+(*Jeq).rows(), (*H).cols());
                m_KKT.middleRows((*H).rows(), (*Jeq).rows()) = (*Jeq);
                m_KKT.conservativeResize((*H).rows()+(*Jeq).rows(), (*H).rows()+(*Jeq).rows());
                
                m_KKTt = m_KKT.transpose();
                m_KKTt.middleRows((*H).rows(), (*Jeq).rows()) = (*Jeq);
                //KKT construction complete
                
                //add the constraint matrix
                m_b.resize((*H).rows()+(*Jeq).rows());
                m_b.segment(0, (*H).rows()) = -g;
                m_b.segment((*H).rows(), (*Jeq).rows()) = (*ceq) - (*Jeq)*x0.head((*H).cols());
                //toMatlab(m_b, "testB.txt");
                
                //Solve and return the newton search direction
                #ifdef GAUSS_PARDISO
                    return m_solver.solve(m_KKTt, m_b);
                #else
                    m_solver.compute(m_KKTt);
            
                    if(m_solver.info()!=Eigen::Success) {
                        // decomposition failed
                        assert(1 == 0);
                        std::cout<<"Decomposition Failed \n";
                        exit(1);
                    }
            
                    return m_solver.solve((*forceVector));
                
                #endif
                
            }
            
        protected:
            
            #ifdef GAUSS_PARDISO
                SolverPardiso<Eigen::SparseMatrix<DataType, Eigen::RowMajor> > m_solver;
            #else
                Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > m_solver;
            #endif

            Eigen::SparseMatrix<DataType, Eigen::RowMajor> m_KKT, m_KKTt;
            Eigen::VectorXx<DataType> m_b;
        };
        
        //Equality constrained newtons method using Gauss
        template <typename Energy, typename Gradient, typename Hessian, typename ConstraintEq,
                  typename JacobianEq, typename Solver, typename PostStepCallback, typename Vector>
        inline bool minimizeBacktracking(Vector &x0, Energy &f, Gradient &g, Hessian &H, ConstraintEq &ceq, JacobianEq &Aeq, Solver &solver, PostStepCallback &pscallback, double tol1 = 1e-5, unsigned int numIterations = UINT_MAX) {
            
            auto linesearch = [&solver, &pscallback](auto &x0, auto &f, auto &g, auto &H, auto &ceq, auto &Aeq, auto  &tol1) {
                return backTrackingLinesearch(x0, f, g, H, ceq, Aeq, solver, pscallback, tol1);
            };
            
            return optimizeWithLineSearch(x0, f, g, H, ceq, Aeq, linesearch, tol1, numIterations);
        }
        
        //functors are annoying but when I have solvers that require initialization they seem to be a necessary evil to avoid reinitilization
        template<typename DataType>
        class NewtonSearchWithBackTracking {
        public:
            
            
            inline NewtonSearchWithBackTracking() { };
            
            template <typename Energy, typename Gradient, typename Hessian, typename ConstraintEq,
            typename JacobianEq, typename PostStepCallback, typename Vector>
            inline bool operator()(Vector &x0, Energy &f, Gradient &g, Hessian &H, ConstraintEq &ceq,
                                   JacobianEq &Aeq, PostStepCallback &pscallback, double tol1=1e-5, unsigned int numIterations= 10000) {
                
                return minimizeBacktracking(x0, f, g, H, ceq, Aeq, m_solver, pscallback, tol1, numIterations);
            }
        protected:
        private:
            
            DirectionNewtonAssembler<DataType> m_solver;
        };
    }
    
}

#endif /* GaussOptimizationAdapters_h */
