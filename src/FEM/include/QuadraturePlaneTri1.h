//
//  QuadratureTri1.h .h
//  Gauss
//
//  Created by David Levin on 9/19/17.
//
//

#ifndef QuadratureTri1_h__h
#define QuadratureTri1_h__h

#include <Quadrature.h>
#include <Energy.h>
#include <UtilitiesMATLAB.h>
#include <ShapeFunctionPlaneLinear.h>

namespace Gauss {
    namespace FEM {
        
        //Integrate energy over 
        template<typename DataType, typename Energy>
        class QuadraturePlaneTri1 : public Energy {
        public:
            
            using Energy::m_qDofs;
            using Energy::m_qDotDofs;
            
            template<typename QDOFList, typename QDotDOFList>
            inline QuadraturePlaneTri1 (Eigen::MatrixXd &V, Eigen::MatrixXi &F,QDOFList &qDOFList, QDotDOFList &qDotDOFList) :
            Energy(V,F,qDOFList, qDotDOFList) { }
            
            inline void getValue(DataType &f, State<DataType> &state) {
                
                DataType w = static_cast<DataType>(Energy::volume());
            
                Eigen::Vector3x<DataType> x;
            
                return w*Energy::getValue(Energy::x(0.5,0.5,0.5).data(), state);
                
            }
            
            template<typename Vector>
            inline void getGradient(Vector &f, const State<DataType> &state) {
                
                DataType w = static_cast<DataType>(Energy::volume());
                
                Eigen::VectorXx<DataType> fInt;
                Eigen::Vector3x<DataType> x;
                fInt.setZero();
                
                
                Energy::getGradient(fInt, Energy::x(0.5,0.5,0.5).data(), state);
                fInt *= w;
                assign(f, fInt, Energy::m_qDofs);
                
            }
            
            template<typename Matrix>
            inline void getHessian(Matrix &H, const State<DataType> &state) {
                DataType w = static_cast<DataType>(Energy::volume());
                
              
                Eigen::MatrixXx<DataType> HInt;
                Eigen::Vector3x<DataType> x;
                HInt.setZero();
                
                
                Energy::getHessian(HInt, Energy::x(0.5,0.5,0.5).data(), state);
                HInt *= w;
                assign(H, HInt, Energy::m_qDofs, Energy::m_qDofs);
                
            }
            
        protected:
            
            
        private:
        };
        
        //Exact integration rule for plane strain triangle mass matrix
        template<typename DataType>
        class QuadratureExact<DataType, EnergyKineticNonLumped<DataType, ShapeFunctionPlaneLinear<DataType> > > :
        public EnergyKineticNonLumped<DataType, ShapeFunctionPlaneLinear<DataType> > {
            
        public:
            
            using EnergyKineticNonLumped<DataType, ShapeFunctionPlaneLinear<DataType> >::m_rho;
            using ShapeFunctionPlaneLinear<DataType>::m_qDotDofs;
            
            template<typename QDOFList, typename QDotDOFList>
            inline QuadratureExact(Eigen::MatrixXd &V, Eigen::MatrixXi &F,QDOFList &qDOFList, QDotDOFList &qDotDOFList) : EnergyKineticNonLumped<DataType, ShapeFunctionPlaneLinear<DataType> >(V,F, qDOFList, qDotDOFList) {
                
                double V0 = EnergyKineticNonLumped<DataType, ShapeFunctionPlaneLinear<DataType> >::volume();
                
                assert(V0 > 0); //tet non inverted in reference config
                
                //setup the mass matrix
                double mass = m_rho*V0;
                double c0 = (1.0/6.0)*mass;
                double c1 = (1.0/12.0)*mass;
                
                //m_massMatrix
                //Point Indices
                unsigned int p0 = 0;
                unsigned int p1 = 3;
                unsigned int p2 = 6;
                
                //Assemble these bad boys //really big 4x4 block matrix
                m_massMatrix.block(p0,p0, 3,3) = c0*Eigen::Matrix<DataType,3,3>::Identity();
                m_massMatrix.block(p0,p1, 3,3) = c1*Eigen::Matrix<DataType,3,3>::Identity();
                m_massMatrix.block(p0,p2, 3,3) = c1*Eigen::Matrix<DataType,3,3>::Identity();
                
                m_massMatrix.block(p1,p0, 3,3) = c1*Eigen::Matrix<DataType,3,3>::Identity();
                m_massMatrix.block(p1,p1, 3,3) = c0*Eigen::Matrix<DataType,3,3>::Identity();
                m_massMatrix.block(p1,p2, 3,3) = c1*Eigen::Matrix<DataType,3,3>::Identity();
                
                m_massMatrix.block(p2,p0, 3,3) = c1*Eigen::Matrix<DataType,3,3>::Identity();
                m_massMatrix.block(p2,p1, 3,3) = c1*Eigen::Matrix<DataType,3,3>::Identity();
                m_massMatrix.block(p2,p2, 3,3) = c0*Eigen::Matrix<DataType,3,3>::Identity();
                
            }
            
            
            inline ~QuadratureExact() { }
            
            //integral rules for things that I want
            template<typename DOFList>
            inline void getValue(DataType &f, const State<DataType> &state) {
                
                
            }
            
            template<typename Vector>
            inline void getGradient(Vector &f, const State<DataType> &state) {
                //do Nothing (I don't think I use this for anything really)
            }
            
            template<typename Matrix>
            inline void getHessian(Matrix &H, const State<DataType> &state) {
                
                assign(H, m_massMatrix, m_qDotDofs, m_qDotDofs);
                
                
            }
            
        protected:
            Eigen::Matrix<DataType, 9,9> m_massMatrix;
            
        private:
        };
        
        //Exact integral for body forces
        template<typename DataType>
        class QuadraturePlaneTri1<DataType, BodyForceGravity<DataType, ShapeFunctionPlaneLinear<DataType> > >  :
            public BodyForceGravity<DataType, ShapeFunctionPlaneLinear<DataType> > {
        public:
            
            using ShapeFunctionPlaneLinear<DataType>::m_qDofs;
            using BodyForceGravity<DataType, ShapeFunctionPlaneLinear<DataType>>::m_rho;
            using BodyForceGravity<DataType, ShapeFunctionPlaneLinear<DataType>>::m_g;
            
                
            template<typename QDOFList, typename QDotDOFList>
            inline QuadraturePlaneTri1 (Eigen::MatrixXd &V, Eigen::MatrixXi &F,QDOFList &qDOFList, QDotDOFList &qDotDOFList) :
            BodyForceGravity<DataType, ShapeFunctionPlaneLinear<DataType> > (V,F,qDOFList, qDotDOFList) { }
            
            inline DataType getValue(State<DataType> &state) {
                
                Eigen::Matrix<double, 9,1> elementF;
                Eigen::Map<Eigen::VectorXd> v0 = mapDOFEigen(*m_qDofs[0], state);
                Eigen::Map<Eigen::VectorXd> v1 = mapDOFEigen(*m_qDofs[1], state);
                Eigen::Map<Eigen::VectorXd> v2 = mapDOFEigen(*m_qDofs[2], state);
                Eigen::Matrix<DataType, 12,1> q;
                q << v0, v1, v2;
                
                getGradient(elementF, state);
                return -q.transpose()*elementF;

            }
            
            template<typename Vector>
            inline void getGradient(Vector &f, const State<DataType> &state) {
                
                Eigen::Vector3d nodalF = m_rho*ShapeFunctionPlaneLinear<DataType>::volume()*(1.0/3.0)*m_g;
                Eigen::Matrix<double, 9,1> elementF = nodalF.replicate(3,1);
                assign(f, elementF, m_qDofs);
                
            }
            
            template<typename Matrix>
            inline void getHessian(Matrix &H, const State<DataType> &state) {
                
            }
            
        protected:
            
            
        private:
        };

    }
}

#endif /* QuadratureTri1_h__h */
