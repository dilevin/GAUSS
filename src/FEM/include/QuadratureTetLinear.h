//
//  QuadratureTetLinear.h
//  Gauss
//
//  Created by David Levin on 12/11/18.
//

#ifndef QuadratureTetLinear_h
#define QuadratureTetLinear_h

namespace Gauss {
    namespace FEM {
        
        template<typename DataType, typename Energy>
        class QuadratureTetLinear : public Energy {
        public:
            
            using Energy::m_qDofs;
            using Energy::m_qDotDofs;
            
            template<typename QDOFList, typename QDotDOFList>
            inline QuadratureTetLinear(Eigen::MatrixXx<DataType> &V, Eigen::MatrixXi &F,QDOFList &qDOFList, QDotDOFList &qDotDOFList) :
            Energy(V,F,qDOFList, qDotDOFList){ }
            
            inline DataType getValue(const State<DataType> &state) {
                DataType w = 0.25*static_cast<DataType>(Energy::volume());
                
                Eigen::Vector3x<DataType> x;
                DataType a,b;
                a = static_cast<DataType>(0.5854101966249685);
                b = static_cast<DataType>(0.1381966011250105);
                
                return w*(Energy::getValue(Energy::x(a,b,b).data(), state) +
                          Energy::getValue(Energy::x(b,a,b).data(), state) +
                          Energy::getValue(Energy::x(b,b,a).data(), state) +
                          Energy::getValue(Energy::x(b,b,b).data(), state));
                
                
            }
            
            template<typename Vector>
            inline void getGradient(Vector &f, const State<DataType> &state) {
                
                DataType w = 0.25*static_cast<DataType>(Energy::volume());
                DataType a = static_cast<DataType>(0.5854101966249685);
                DataType b = static_cast<DataType>(0.1381966011250105);
                
                Eigen::VectorXx<DataType> fInt, fPoint;
                Eigen::Vector3x<DataType> x;
                fInt.setZero();
                
                Energy::getGradient(fPoint, Energy::x(a,b, b).data(), state);
                fInt = fPoint*w;
                
                Energy::getGradient(fPoint, Energy::x(b,a, b).data(), state);
                fInt += fPoint*w;
                
                Energy::getGradient(fPoint, Energy::x(b,b, a).data(), state);
                fInt += fPoint*w;
                
                Energy::getGradient(fPoint, Energy::x(b,b, b).data(), state);
                fInt += fPoint*w;
                
                assign(f, fInt, Energy::m_qDofs);
                
                
            }
            
            template<typename Matrix>
            inline void getHessian(Matrix &H, const State<DataType> &state) {
                DataType w = 0.25*static_cast<DataType>(Energy::volume());
                DataType a = static_cast<DataType>(0.5854101966249685);
                DataType b = static_cast<DataType>(0.1381966011250105);
                
                
                Eigen::MatrixXx<DataType> HInt, HPoint;
                Eigen::Vector3x<DataType> x;
                HInt.setZero();
                
                Energy::getHessian(HPoint, Energy::x(a, b,b).data(), state);
                HInt = w*HPoint;
                
                Energy::getHessian(HPoint, Energy::x(b, a,b).data(), state);
                HInt += w*HPoint;
                
                Energy::getHessian(HPoint, Energy::x(b, b,a).data(), state);
                HInt += w*HPoint;
                
                Energy::getHessian(HPoint, Energy::x(b, b,b).data(), state);
                HInt += w*HPoint;
                
                assign(H, HInt, Energy::m_qDofs, Energy::m_qDofs);
                
            }
            
        protected:
        private:
        };
    }
}

#endif /* QuadratureTetLinear_h */
