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
        class QuadratureHex8 : public Energy {
        public:
            
            using Energy::m_dx;
            using Energy::m_x0;
            using Energy::m_qDofs;
            using Energy::m_qDotDofs;
            
            template<typename QDOFList, typename QDotDOFList>
            inline QuadratureHex8(Eigen::MatrixXd &V, Eigen::MatrixXi &F,QDOFList &qDOFList, QDotDOFList &qDotDOFList) :
            Energy(V,F,qDOFList, qDotDOFList){ }
            
            inline double getValue(State<DataType> &state) {
                DataType w = static_cast<DataType>(Energy::volume())/8.0;
                DataType alpha = static_cast<DataType>(sqrt(1.0/3.0));
                
                double energy = 0.0;
                Eigen::Vector3x<DataType> x;
                
                energy += w*Energy::getValue(Energy::x(-alpha,-alpha,-alpha).data(), state);
                
                energy += w*Energy::getValue(Energy::x(alpha,-alpha,-alpha).data(), state);
                
                energy += w*Energy::getValue(Energy::x(alpha,alpha,-alpha).data(), state);
                
                energy += w*Energy::getValue(Energy::x(-alpha,alpha,-alpha).data(), state);
                
                energy += w*Energy::getValue(Energy::x(-alpha,-alpha,alpha).data(), state);
                
                energy += w*Energy::getValue(Energy::x(alpha,-alpha,alpha).data(), state);
                
                energy += w*Energy::getValue(Energy::x(alpha,alpha,alpha).data(), state);
                
                energy += w*Energy::getValue(Energy::x(-alpha,alpha,alpha).data(), state);
                
                return energy; 
            }
            
            template<typename Vector>
            inline void getGradient(Vector &f, const State<DataType> &state) {
            
                DataType w = static_cast<DataType>(Energy::volume())/8.0;
                DataType alpha = static_cast<DataType>(sqrt(1.0/3.0));
            
                Eigen::VectorXx<DataType> fInt;
                Eigen::Vector3x<DataType> x;
                fInt.setZero();
                
                
                Energy::getGradient(fInt, Energy::x(-alpha,-alpha,-alpha).data(), state);
                fInt *= w;
                assign(f, fInt, Energy::m_qDofs);
                
                Energy::getGradient(fInt, Energy::x(alpha,-alpha,-alpha).data(), state);
                fInt *= w;
                assign(f, fInt, Energy::m_qDofs);
                
                Energy::getGradient(fInt, Energy::x(alpha,alpha,-alpha).data(), state);
                fInt *= w;
                assign(f, fInt, Energy::m_qDofs);
                
                Energy::getGradient(fInt, Energy::x(-alpha,alpha,-alpha).data(), state);
                fInt *= w;
                assign(f, fInt, Energy::m_qDofs);
                
                Energy::getGradient(fInt, Energy::x(-alpha,-alpha,alpha).data(), state);
                fInt *= w;
                assign(f, fInt, Energy::m_qDofs);
                
                Energy::getGradient(fInt, Energy::x(alpha,-alpha,alpha).data(), state);
                fInt *= w;
                assign(f, fInt, Energy::m_qDofs);
                
                Energy::getGradient(fInt, Energy::x(alpha,alpha,alpha).data(), state);
                fInt *= w;
                assign(f, fInt, Energy::m_qDofs);
                
                Energy::getGradient(fInt, Energy::x(-alpha,alpha,alpha).data(), state);
                fInt *= w;
                assign(f, fInt, Energy::m_qDofs);
                
            }
            
            template<typename Matrix>
            inline void getHessian(Matrix &H, const State<DataType> &state) {
                DataType w = static_cast<DataType>(Energy::volume())/8.0;

                DataType alpha = static_cast<DataType>(sqrt(1.0/3.0));
                
                Eigen::MatrixXx<DataType> HInt;
                Eigen::Vector3x<DataType> x;
                HInt.setZero();
                
                
                Energy::getHessian(HInt, Energy::x(-alpha,-alpha,-alpha).data(), state);
                HInt *= w;
                assign(H, HInt, Energy::m_qDofs, Energy::m_qDofs);
                
                Energy::getHessian(HInt, Energy::x(alpha,-alpha,-alpha).data(), state);
                HInt *= w;
                assign(H, HInt, Energy::m_qDofs, Energy::m_qDofs);
                
                Energy::getHessian(HInt, Energy::x(alpha,alpha,-alpha).data(), state);
                HInt *= w;
                assign(H, HInt, Energy::m_qDofs, Energy::m_qDofs);
                
                Energy::getHessian(HInt, Energy::x(-alpha,alpha,-alpha).data(), state);
                HInt *= w;
                assign(H, HInt, Energy::m_qDofs, Energy::m_qDofs);
                
                Energy::getHessian(HInt, Energy::x(-alpha,-alpha,alpha).data(), state);
                HInt *= w;
                assign(H, HInt, Energy::m_qDofs, Energy::m_qDofs);
                
                Energy::getHessian(HInt, Energy::x(alpha,-alpha,alpha).data(), state);
                HInt *= w;
                assign(H, HInt, Energy::m_qDofs, Energy::m_qDofs);
                
                Energy::getHessian(HInt, Energy::x(alpha,alpha,alpha).data(), state);
                HInt *= w;
                assign(H, HInt, Energy::m_qDofs, Energy::m_qDofs);
                
                Energy::getHessian(HInt, Energy::x(-alpha,alpha,alpha).data(), state);
                HInt *= w;
                assign(H, HInt, Energy::m_qDofs, Energy::m_qDofs);
            }
            
        protected:
        private:
        };
    }
}

#endif /* QuadratureHex8_h */
