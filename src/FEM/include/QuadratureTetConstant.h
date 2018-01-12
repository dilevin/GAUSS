//
//  QuadratureTetConstant.h
//  Gauss
//
//  Created by David Levin on 11/7/17.
//
//

#ifndef QuadratureTetConstant_h
#define QuadratureTetConstant_h

namespace Gauss {
    namespace FEM {
        
        template<typename DataType, typename Energy>
        class QuadratureTetConstant : public Energy {
        public:
            
            using Energy::m_qDofs;
            using Energy::m_qDotDofs;
            
            template<typename QDOFList, typename QDotDOFList>
            inline QuadratureTetConstant(Eigen::MatrixXd &V, Eigen::MatrixXi &F,QDOFList &qDOFList, QDotDOFList &qDotDOFList) :
            Energy(V,F,qDOFList, qDotDOFList){ }
            
            inline DataType getValue(const State<DataType> &state) {
                DataType w = static_cast<DataType>(Energy::volume());
                
                Eigen::Vector3x<DataType> x;
                
                return w*Energy::getValue(Energy::x(0.0,0.0,0.0).data(), state);
                
            }
            
            template<typename Vector>
            inline void getGradient(Vector &f, const State<DataType> &state) {
                
                DataType w = static_cast<DataType>(Energy::volume());
                
                Eigen::VectorXx<DataType> fInt;
                Eigen::Vector3x<DataType> x;
                fInt.setZero();
                
                Energy::getGradient(fInt, Energy::x(0.0,0.0,0.0).data(), state);
                fInt *= w;
                assign(f, fInt, Energy::m_qDofs);
                
                
            }
            
            template<typename Matrix>
            inline void getHessian(Matrix &H, const State<DataType> &state) {
                DataType w = static_cast<DataType>(Energy::volume());
                
                Eigen::MatrixXx<DataType> HInt;
                Eigen::Vector3x<DataType> x;
                HInt.setZero();
                
                Energy::getHessian(HInt, Energy::x(0.0, 0.0,0.0).data(), state);
                HInt *= w;
                assign(H, HInt, Energy::m_qDofs, Energy::m_qDofs);
                
            }
            
        protected:
        private:
        };
    }
}

#endif /* QuadratureTetConstant_h */
