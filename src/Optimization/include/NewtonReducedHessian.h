//
//  NewtonReducedHessian.h
//  Gauss
//
//  Created by David Levin on 6/25/18.
//

#ifndef NewtonReducedHessian_h
#define NewtonReducedHessian_h

#include <Newton.h>
#include <GaussOptimizationAdapters.h>
#include <igl/orth.h>

namespace Gauss {
    namespace Optimization {
        //Compute search direction using a reduced form of the stiffness matrix
        template<typename DataType>
        class DirectionNewtonAssemblerReduced : public DirectionNewtonAssembler<DataType> {
            
        public:
            inline DirectionNewtonAssemblerReduced() : DirectionNewtonAssembler<DataType>()
            {
                m_numBasis = 0;
                m_maxBasis = 100;
                m_first = true;
            }
            
            void setU(Eigen::MatrixXx<DataType> &U) { m_U = U; }
            
            //takes in assembled matrices and returns the newton step direction
            template<typename Hessian, typename Gradient, typename Ceq, typename JacobianEq, typename Vector>
            inline auto & operator()(Hessian &H, Gradient &g, Ceq &ceq, JacobianEq &Jeq, Vector &x0) {
                
                if(m_first) {
                    m_first = false;
                    /*Eigen::SparseMatrix<DataType, Eigen::RowMajor> Ident = (*H);
                    Ident.setIdentity();
                    auto eigs = generalizedEigenvalueProblem((*H), (*H), m_maxBasis,   1e-6);
                    m_U = eigs.first;
                    m_first = false;*/
                
                    
                    /*m_g0 = g;
                    m_U.conservativeResize((*H).rows(), m_maxBasis);
                    m_p = DirectionNewtonAssembler<DataType>::operator()(H, g, ceq, Jeq, x0);
                    m_U.rightCols(m_maxBasis-2) = m_U.leftCols(m_maxBasis - 2);
                    m_U.col(0) = m_p.head(m_U.rows());
                    m_U.col(1) = g;
                    m_numBasis += 2;*/
                    m_p = DirectionNewtonAssembler<DataType>::operator()(H, g, ceq, Jeq, x0);
                    return m_p;
                }
                
                //I think this is all wrong
                //step 1 solve unconstrained
                m_H = m_U.transpose()*(*H)*m_U;
                m_p.head(m_U.rows()) = m_U*m_H.ldlt().solve(m_U.transpose()*g);
                
                //step 2 solve for lagrange multipliers
                m_J = (*Jeq)*m_U;
                m_L = m_J*m_H.ldlt().solve(m_J.transpose());
                //m_lambda = m_L.ldlt().solve(m_J*m_U.transpose()*m_p.head(m_U.rows()));
                
                //step 3 build constrained problem
                m_p.head(m_U.rows()) = -m_p.head(m_U.rows()) ;// + m_U*m_H.ldlt().solve(m_J.transpose()*m_lambda);
                m_p.tail((*Jeq).rows()) = m_lambda;
                
                //update basis
                /*auto temp = m_U.leftCols(m_maxBasis - 1);
                m_U.rightCols(m_maxBasis-1) = temp;
                m_U.col(0) = g;
                m_numBasis = std::min(m_maxBasis, ++m_numBasis);
                
                Eigen::MatrixXx<DataType> Q;
                igl::orth(m_U.leftCols(m_numBasis), Q);
                m_numBasis = Q.cols();
                m_U.leftCols(m_numBasis) = Q;
                
                //std::cout<<m_U.leftCols(m_numBasis)<<"\n";
                //solve for m_p
                //step 1 solve unconstrained
                m_H = (m_U.leftCols(m_numBasis).transpose()*(*H)*m_U.leftCols(m_numBasis));
                m_p.head(m_U.rows()) = m_U.leftCols(m_numBasis)*m_H.ldlt().solve(m_U.leftCols(m_numBasis).transpose()*g);
                
                //step 2 solve for lagrange multipliers
                m_J = (*Jeq)*m_U.leftCols(m_numBasis);
                m_L = m_J*m_H.ldlt().solve(m_J.transpose());
                m_lambda = m_L.ldlt().solve(m_J*m_U.leftCols(m_numBasis).transpose()*m_p.head(m_U.rows()));
                
                //step 3 build constrained problem
                m_p.head(m_U.rows()) = -m_p.head(m_U.rows());// + m_U.leftCols(m_numBasis)*m_H.ldlt().solve(m_J.transpose()*m_lambda);
                m_p.tail((*Jeq).rows()) = m_lambda;
                temp = m_U.leftCols(m_maxBasis - 1);
                m_U.rightCols(m_maxBasis-1) = temp;
                m_U.col(0) = m_p.head(m_U.rows());
                m_numBasis = std::min(m_maxBasis, ++m_numBasis);
                m_g0 = g;*/
                return m_p;
            }
            
        protected:
            
            unsigned int m_numBasis;
            unsigned int m_maxBasis;
            bool m_first;

            Eigen::SimplicialLDLT<Eigen::MatrixXx<DataType> > m_denseSolver;
            Eigen::MatrixXx<DataType> m_H, m_L, m_U, m_J;
            Eigen::VectorXx<DataType> m_g0, m_p, m_lambda;
            
        };

        //Newton's Solve using reduced basis
        template<typename DataType>
        class NewtonSearchWithBackTrackingReduced {
        public:
            
            
            inline NewtonSearchWithBackTrackingReduced() { };
            
            template <typename Energy, typename Gradient, typename Hessian, typename ConstraintEq,
            typename JacobianEq, typename PostStepCallback, typename Vector>
            inline bool operator()(Vector &x0, Energy &f, Gradient &g, Hessian &H, ConstraintEq &ceq,
                                   JacobianEq &Aeq, PostStepCallback &pscallback, double tol1=1e-5, unsigned int numIterations= 10000) {
                
                return minimizeBacktracking(x0, f, g, H, ceq, Aeq, m_solver, pscallback, tol1, numIterations);
            }
            
            void setU(Eigen::MatrixXx<DataType> U) { m_solver.setU(U); }
        protected:
        private:
            
            DirectionNewtonAssemblerReduced<DataType> m_solver;
        };
    }
}
#endif /* NewtonReducedHessian_h */
