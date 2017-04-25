//
//  Quadrature.h
//  Gauss
//
//  Created by David Levin on 3/13/17.
//
//

#ifndef Quadrature_h
#define Quadrature_h

#include <State.h>

namespace Gauss {
    namespace FEM {
        template<typename DataType, typename Energy>
        class QuadratureNone : public Energy {
        public:
            
            inline void getValue(DataType &f, State<DataType> &state) { }

            template<typename Vector, typename DOFList>
            inline void getGradient(Vector &f, State<DataType> &state) { }
            
            template<typename Matrix, typename DOFList>
            inline void getHessian(Matrix &H, State<DataType> &state, const DOFList &q) { }

        };
        
        //Specialize for exact integrals
        template<typename DataType, typename Energy>
        class QuadratureExact : public Energy {
        public:
            template<typename DOFList>
            QuadratureExact(Eigen::MatrixXd &V, Eigen::MatrixXi &F, DOFList &dofList) : Energy(V,F, dofList) {
                std::cout<<"Do nothing integration rule \n";
            }
            ~QuadratureExact() { }
        protected:
        private:
        };
        
    }
}



#endif /* Quadrature_h */
