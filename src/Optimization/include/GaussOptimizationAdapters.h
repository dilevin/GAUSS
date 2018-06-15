//
//  GaussOptimizationAdapters.h
//  Gauss
//
//  Created by David Levin on 6/13/18.
//

#ifndef GaussOptimizationAdapters_h
#define GaussOptimizationAdapters_h

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
        
            template<typename Vector, template Func>
            inline auto operator()(Vector &qDot) {
                
                //init assembler
                ASSEMBLEMATINIT(m_assembler, m_m, m_n);
                return this->operator()(qDot, m_assemler);
                ASSEMBLEEND(m_assembler);
                
            }
            
        protected:
            Assembler<DataType> &m_assembler;
            Func &m_func;
            unsigned int m_m, m_n; //assembler sizes
            
        };

        //Newtons Step Direction using Pardiso
        template<typename DataType>
        class DirectionNewtonPardiso {
          
        public:
            inline DirectionNewtonPardiso() {
                
            }
            
            template<typename HessianFunc, typename GradientFunc, typename Ceq, typename JacobianEq>
            inline auto operator()(HessianFunc &H, GradientFunc &g, Ceq &cdq, JacobianEq &Jeq) {
                
                //compute newton search direction using proper hessian
                
                //return newton search direction
            }
        protected:
            
            SolverPardiso<Eigen::SparseMatrix<DataType, Eigen::RowMajor> > m_pardiso;
            AssemblerParallel<DataType, AssemblerSparseMatrixEigen<DataType> > m_mat;
            AssemblerParallel<DataType, AssemblerVectorEigen<DataType> > m_vec;
        };
    }
}

#endif /* GaussOptimizationAdapters_h */
