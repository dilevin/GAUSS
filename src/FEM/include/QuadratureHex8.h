//
//  QuadratureHex8.h
//  Gauss
//
//  Created by David Levin on 5/9/17.
//
//

#ifndef QuadratureHex8_h
#define QuadratureHex8_h

namespace Gauss {
    namespace FEM {
        
        template<typename DataType, typename Energy>
        class QuadratureHex8 : Energy {
        public:
            
            using Energy::m_dx;
            using Energy::m_x0;
            
            template<typename QDOFList, typename QDotDOFList>
            inline QuadratureHex8(Eigen::MatrixXd &V, Eigen::MatrixXi &F,QDOFList &qDOFList, QDotDOFList &qDotDOFList) :
            Energy(V,F,qDOFList, qDotDOFList){ }
            
            inline void getValue(DataType &f, State<DataType> &state) {
                
            }
            
            template<typename Vector, typename DOFList>
            inline void getGradient(Vector &f, State<DataType> &state) {
            
                DataType w = static_cast<DataType>(1.0);
                DataType alpha = static_cast<DataType>(sqrt(1.0/3.0));
            
                Eigen::VectorXx<DataType> fInt;
                Eigen::Vector3x<DataType> x;
                fInt.setZero();
                
                
                Energy::getGradient(fInt, Energy::x(-alpha,-alpha,-alpha).data(), state);
                fInt *= w;
                assign(f, fInt, Energy::m_qDof);
                
                Energy::getGradient(fInt, Energy::x(alpha,-alpha,-alpha).data(), state);
                fInt *= w;
                assign(f, fInt, Energy::m_qDof);
                
                Energy::getGradient(fInt, Energy::x(alpha,alpha,-alpha).data(), state);
                fInt *= w;
                assign(f, fInt, Energy::m_qDof);
                
                Energy::getGradient(fInt, Energy::x(alpha,-alpha,-alpha).data(), state);
                fInt *= w;
                assign(f, fInt, Energy::m_qDof);
                
                Energy::getGradient(fInt, Energy::x(-alpha,-alpha,alpha).data(), state);
                fInt *= w;
                assign(f, fInt, Energy::m_qDof);
                
                Energy::getGradient(fInt, Energy::x(alpha,-alpha,alpha).data(), state);
                fInt *= w;
                assign(f, fInt, Energy::m_qDof);
                
                Energy::getGradient(fInt, Energy::x(alpha,alpha,alpha).data(), state);
                fInt *= w;
                assign(f, fInt, Energy::m_qDof);
                
                Energy::getGradient(fInt, Energy::x(alpha,-alpha,alpha).data(), state);
                fInt *= w;
                assign(f, fInt, Energy::m_qDof);
                
                
                
            }
            
            template<typename Matrix, typename DOFList>
            inline void getHessian(Matrix &H, State<DataType> &state, const DOFList &q) {
                DataType w = static_cast<DataType>(1.0);
                DataType alpha = static_cast<DataType>(sqrt(1.0/3.0));
                
                Eigen::MatrixXx<DataType> HInt;
                Eigen::Vector3x<DataType> x;
                HInt.setZero();
                
                
                Energy::getHessian(HInt, Energy::x(-alpha,-alpha,-alpha).data(), state);
                HInt *= w;
                assign(H, HInt, Energy::m_qDof, Energy::m_qDof);
                
                Energy::getHessian(HInt, Energy::x(alpha,-alpha,-alpha).data(), state);
                HInt *= w;
                assign(H, HInt, Energy::m_qDof, Energy::m_qDof);
                
                Energy::getHessian(HInt, Energy::x(alpha,alpha,-alpha).data(), state);
                HInt *= w;
                assign(H, HInt, Energy::m_qDof, Energy::m_qDof);
                
                Energy::getHessian(HInt, Energy::x(alpha,-alpha,-alpha).data(), state);
                HInt *= w;
                assign(H, HInt, Energy::m_qDof, Energy::m_qDof);
                
                Energy::getHessian(HInt, Energy::x(-alpha,-alpha,alpha).data(), state);
                HInt *= w;
                assign(H, HInt, Energy::m_qDof, Energy::m_qDof);
                
                Energy::getHessian(HInt, Energy::x(alpha,-alpha,alpha).data(), state);
                HInt *= w;
                assign(H, HInt, Energy::m_qDof, Energy::m_qDof);
                
                Energy::getHessian(HInt, Energy::x(alpha,alpha,alpha).data(), state);
                HInt *= w;
                assign(H, HInt, Energy::m_qDof, Energy::m_qDof);
                
                Energy::getHessian(HInt, Energy::x(alpha,-alpha,alpha).data(), state);
                HInt *= w;
                assign(H, HInt, Energy::m_qDof, Energy::m_qDof);
            }
            
        protected:
        private:
        };
    }
}

#endif /* QuadratureHex8_h */
