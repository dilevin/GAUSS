//
//  Element.h
//  Gauss
//
//  Created by David Levin on 3/13/17.
//
//

#ifndef Element_h
#define Element_h

#include <Quadrature.h>
#include <QuadratureExact.h>
#include <ShapeFunction.h>
#include <Energy.h>

//Elements allow the evaluation of functions over the domain.
//maybe need to do this with trait classes
namespace Gauss {
    namespace FEM {
        template<   typename DataType, unsigned int N,
                    template<typename Type, typename Energy> class QuadratureRuleT,
                    template<typename Type, typename Energy> class QuadratureRuleU,
                    template<typename Type, typename Func> class KineticEnergy,
                    template<typename Type, typename Func> class PotentialEnergy,
                    template<typename Type, typename Func> class BodyForce,
                    template<typename Type> class ShapeFunction >
        class ElementBase  : public QuadratureRuleU<DataType, PotentialEnergy<DataType, ShapeFunction<DataType> > >,
                         public QuadratureRuleT<DataType, KineticEnergy<DataType, ShapeFunction<DataType> > >,
                         public QuadratureRuleU<DataType, BodyForce<DataType, ShapeFunction<DataType> > > {
        public:
                             
            using QuadratureT = QuadratureRuleT<DataType, KineticEnergy<DataType, ShapeFunction<DataType> > >;
            using QuadratureU = QuadratureRuleU<DataType, PotentialEnergy<DataType, ShapeFunction<DataType> > >;
            using QuadratureBF = QuadratureRuleU<DataType, BodyForce<DataType, ShapeFunction<DataType> > >;

            template<typename QDOFList, typename QDotDOFList>
            ElementBase(Eigen::MatrixXd &V, Eigen::MatrixXi &F, QDOFList qDOFList, QDotDOFList qDotDOFList)  : QuadratureU(V,F, qDOFList, qDotDOFList), QuadratureT(V,F, qDOFList, qDotDOFList), QuadratureBF(V,F, qDOFList, qDotDOFList) {
            
            }
            ~ElementBase() { }
            
            static constexpr unsigned int numDOFs() { return N; }
                             
            void getKineticEnergy(DataType &scalar, State<DataType> &state) {
                QuadratureT::getValue(scalar, state);
            }
            
                             
            void getPotentialEnergy(DataType &scalar, State<DataType> &state) {
                QuadratureU::getValue(scalar, state);
            }
            
            template<typename Vector>
            void getForce(Vector &f, const State<DataType> &state) {
                QuadratureU::getGradient(f, state);
                getBodyForce(f, state);
            }
            
            template<typename Matrix>
            void getMassMatrix(Matrix &M, const State<DataType> &state) {
                QuadratureT::getHessian(M, state);
            }
                             
            template<typename Matrix>
            void getStiffnessMatrix(Matrix &H, const State<DataType> &state) {
                  QuadratureU::getHessian(H, state);
            }
            
            template<typename Vector>
            inline void getBodyForce(Vector &f, const State<DataType> &state) {
                 //Integrate Body Force (send it density function)
                QuadratureBF::setParams(QuadratureT::getDensity(), 0.0, -9.8, 0.0); //temporary, make more efficient later
                QuadratureBF::getGradient(f, state);
            }
                             
        protected:
                             
                        
        private:
                             
        };
            
        
        //convenience definition for shared potential and kinetic energy quadrature rules
        template<   typename DataType, unsigned int N,
        template<typename Type, typename Energy> class QuadratureRule,
        template<typename Type, typename Func> class KineticEnergy,
        template<typename Type, typename Func> class PotentialEnergy,
        template<typename Type, typename Func> class BodyForce,
        template<typename Type> class ShapeFunction >
        using Element = ElementBase<DataType, N, QuadratureRule, QuadratureRule, KineticEnergy, PotentialEnergy, BodyForce, ShapeFunction>;
        
    }
}
#endif /* Element_h */
