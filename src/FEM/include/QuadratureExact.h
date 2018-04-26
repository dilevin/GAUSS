//
//  QuadratureExact.h
//  Gauss
//
//  Created by David Levin on 3/14/17.
//
//

#ifndef QuadratureExact_h
#define QuadratureExact_h

#include <Quadrature.h>
#include <Energy.h>
#include <UtilitiesMATLAB.h>

//Some exact quadrature rules
namespace Gauss {
    namespace FEM {
        //For linear tetrahedra
        template<typename DataType, typename ShapeFunction>
        class QuadratureExact<DataType, EnergyKineticNonLumped<DataType, ShapeFunction> > :
                                    public EnergyKineticNonLumped<DataType, ShapeFunction> {
                    
        public:
            
            using EnergyKineticNonLumped<DataType, ShapeFunction>::m_rho;
            using ShapeFunction::m_qDotDofs;
                                        
            template<typename QDOFList, typename QDotDOFList>
            inline QuadratureExact(Eigen::MatrixXd &V, Eigen::MatrixXi &F,QDOFList &qDOFList, QDotDOFList &qDotDOFList) : EnergyKineticNonLumped<DataType, ShapeFunction>(V,F, qDOFList, qDotDOFList) {
                
                //need the volume of this tetrahedron
                Eigen::Matrix<DataType, 3,3> V0;
                V0 << Vert(1,0) - Vert(0,0), Vert(2,0) - Vert(0,0), Vert(3,0) - Vert(0 ,0),
                Vert(1,1) - Vert(0,1), Vert(2,1) - Vert(0,1), Vert(3,1) - Vert(0,1),
                Vert(1,2) - Vert(0,2), Vert(2,2) - Vert(0,2), Vert(3,2) - Vert(0,2);
 
                m_V0 = (1.0/6.0)*V0.determinant();
                
                if(m_V0 <= 0) {
                    std::cout<<"Inverted element detected \n";
                    exit(1);
                }
                assert(m_V0 > 0); //tet non inverted in reference config
                
                //setup the mass matrix
                double mass = m_rho*m_V0;
                double c0 = (1.0/10.0)*mass;
                double c1 = (1.0/20.0)*mass;
                
                //m_massMatrix
                //Point Indices
                unsigned int p0 = 0;
                unsigned int p1 = 3;
                unsigned int p2 = 6;
                unsigned int p3 = 9;
                
                //Assemble these bad boys //really big 4x4 block matrix
                m_massMatrix.block(p0,p0, 3,3) = c0*Eigen::Matrix<DataType,3,3>::Identity();
                m_massMatrix.block(p0,p1, 3,3) = c1*Eigen::Matrix<DataType,3,3>::Identity();
                m_massMatrix.block(p0,p2, 3,3) = c1*Eigen::Matrix<DataType,3,3>::Identity();
                m_massMatrix.block(p0,p3, 3,3) = c1*Eigen::Matrix<DataType,3,3>::Identity();
                
                m_massMatrix.block(p1,p0, 3,3) = c1*Eigen::Matrix<DataType,3,3>::Identity();
                m_massMatrix.block(p1,p1, 3,3) = c0*Eigen::Matrix<DataType,3,3>::Identity();
                m_massMatrix.block(p1,p2, 3,3) = c1*Eigen::Matrix<DataType,3,3>::Identity();
                m_massMatrix.block(p1,p3, 3,3) = c1*Eigen::Matrix<DataType,3,3>::Identity();
                
                m_massMatrix.block(p2,p0, 3,3) = c1*Eigen::Matrix<DataType,3,3>::Identity();
                m_massMatrix.block(p2,p1, 3,3) = c1*Eigen::Matrix<DataType,3,3>::Identity();
                m_massMatrix.block(p2,p2, 3,3) = c0*Eigen::Matrix<DataType,3,3>::Identity();
                m_massMatrix.block(p2,p3, 3,3) = c1*Eigen::Matrix<DataType,3,3>::Identity();
                
                m_massMatrix.block(p3,p0, 3,3) = c1*Eigen::Matrix<DataType,3,3>::Identity();
                m_massMatrix.block(p3,p1, 3,3) = c1*Eigen::Matrix<DataType,3,3>::Identity();
                m_massMatrix.block(p3,p2, 3,3) = c1*Eigen::Matrix<DataType,3,3>::Identity();
                m_massMatrix.block(p3,p3, 3,3) = c0*Eigen::Matrix<DataType,3,3>::Identity();
                
            }
                                        
                                        
            inline ~QuadratureExact() { }
            
            //integral rules for things that I wan
            inline DataType getValue(const State<DataType> &state) {
            
                Eigen::Map<Eigen::VectorXd> v0 = mapDOFEigen(*m_qDotDofs[0], state);
                Eigen::Map<Eigen::VectorXd> v1 = mapDOFEigen(*m_qDotDofs[1], state);
                Eigen::Map<Eigen::VectorXd> v2 = mapDOFEigen(*m_qDotDofs[2], state);
                Eigen::Map<Eigen::VectorXd> v3 = mapDOFEigen(*m_qDotDofs[3], state);
                Eigen::Matrix<DataType, 12,1> qDot;
                qDot << v0, v1, v2, v3;
                
                return 0.5*qDot.transpose()*m_massMatrix*qDot;
            }
            
            template<typename Vector>
            inline void getGradient(Vector &f, const State<DataType> &state) {
                //do Nothing (I don't think I use this for anything really)
            }
            
            template<typename Matrix>
            inline void getHessian(Matrix &H, const State<DataType> &state) {
            
                assign(H, m_massMatrix, m_qDotDofs, m_qDotDofs);

                
                //Only need two numbers for element mass matrix
                /*double mass = m_rho*m_V0;
                double c0 = (1.0/10.0)*mass;
                double c1 = (1.0/20.0)*mass;
                int status = 1;
                
                //Point Indices
                unsigned int p0 = 3*m_pointIndices.get(0);
                unsigned int p1 = 3*m_pointIndices.get(1);
                unsigned int p2 = 3*m_pointIndices.get(2);
                unsigned int p3 = 3*m_pointIndices.get(3);
                
                //Assemble these bad boys //really big 4x4 block matrix
                status *= LinearAlgebra::assemble(matrix, c0, 3, offset+p0, offset+p0, false);
                status *= LinearAlgebra::assemble(matrix, c1, 3, offset+p0, offset+p1, false);
                status *= LinearAlgebra::assemble(matrix, c1, 3, offset+p0, offset+p2, false);
                status *= LinearAlgebra::assemble(matrix, c1, 3, offset+p0, offset+p3, false);
                
                status *= LinearAlgebra::assemble(matrix, c1, 3, offset+p1, offset+p0, false);
                status *= LinearAlgebra::assemble(matrix, c0, 3, offset+p1, offset+p1, false);
                status *= LinearAlgebra::assemble(matrix, c1, 3, offset+p1, offset+p2, false);
                status *= LinearAlgebra::assemble(matrix, c1, 3, offset+p1, offset+p3, false);
                
                status *= LinearAlgebra::assemble(matrix, c1, 3, offset+p2, offset+p0, false);
                status *= LinearAlgebra::assemble(matrix, c1, 3, offset+p2, offset+p1, false);
                status *= LinearAlgebra::assemble(matrix, c0, 3, offset+p2, offset+p2, false);
                status *= LinearAlgebra::assemble(matrix, c1, 3, offset+p2, offset+p3, false);
                
                status *= LinearAlgebra::assemble(matrix, c1, 3, offset+p3, offset+p0, false);
                status *= LinearAlgebra::assemble(matrix, c1, 3, offset+p3, offset+p1, false);
                status *= LinearAlgebra::assemble(matrix, c1, 3, offset+p3, offset+p2, false);
                status *= LinearAlgebra::assemble(matrix, c0, 3, offset+p3, offset+p3, false);*/

            }
            
        protected:
            Eigen::Matrix<DataType, 12,12> m_massMatrix;
            double m_V0; //volume
            
        private:
        };
        
        template<typename DataType, typename ShapeFunction>
        class QuadratureExact<DataType, EnergyLinearElasticity<DataType, ShapeFunction> > :
        public EnergyLinearElasticity<DataType, ShapeFunction> {
            
        public:
            
            using ShapeFunction::dphi;
            
            using EnergyLinearElasticity<DataType, ShapeFunction>::m_E;
            using EnergyLinearElasticity<DataType, ShapeFunction>::m_mu;
            using ShapeFunction::m_qDofs;
            using ShapeFunction::m_qDotDofs;
            
            template<typename QDOFList, typename QDotDOFList>
            inline QuadratureExact(Eigen::MatrixXd &V, Eigen::MatrixXi &F, QDOFList &qDOFList, QDotDOFList &qDotDOFList) : EnergyLinearElasticity<DataType, ShapeFunction>(V,F, qDOFList, qDotDOFList) {
                
                //need the volume of this tetrahedron
                Eigen::Matrix<DataType, 3,3> V0;
                V0 << Vert(1,0) - Vert(0,0), Vert(2,0) - Vert(0,0), Vert(3,0) - Vert(0,0),
                Vert(1,1) - Vert(0,1), Vert(2,1) - Vert(0,1), Vert(3,1) - Vert(0,1),
                Vert(1,2) - Vert(0,2), Vert(2,2) - Vert(0,2), Vert(3,2) - Vert(0,2);
                
                m_V0 = (1.0/6.0)*V0.determinant();
                
                assert(m_V0 > 0); //chewck tet is not inverted in reference config
                
                //Compute B and C matrices then form per-element K
                //Form of B if we assume DOFS are tiles as [x,y,z, x,y,z ....]'
                //
                Eigen::Matrix<DataType, 6, 6> C;
                C.setZero();
                C(0,0) = 1.0-m_mu;
                C(0,1) = m_mu;
                C(0,2) = m_mu;
                C(1,0) = m_mu;
                C(1,1) = 1.0-m_mu;
                C(1,2) = m_mu;
                C(2,0) = m_mu;
                C(2,1) = m_mu;
                C(2,2) = 1.0-m_mu;
                C(3,3) = 0.5*(1.0-2.0*m_mu);
                C(4,4) = 0.5*(1.0-2.0*m_mu);
                C(5,5) = 0.5*(1.0-2.0*m_mu);
                C *= (m_E/((1.0+m_mu)*(1.0-2.0*m_mu)));
                
                std::array<DataType, 3> x = {{0,0,0}};
                std::array<DataType, 3> dphi0 = this->template dphi<0>(x.data());
                std::array<DataType, 3> dphi1 = this->template dphi<1>(x.data());
                std::array<DataType, 3> dphi2 = this->template dphi<2>(x.data());
                std::array<DataType, 3> dphi3 = this->template dphi<3>(x.data());
                
                Eigen::Matrix<DataType, 6,12> B;
                B <<    dphi0[0],  0,  0,  dphi1[0],   0,  0,  dphi2[0],   0,  0,  dphi3[0],   0,   0,
                        0,      dphi0[1],  0,  0,  dphi1[1],   0,  0,  dphi2[1],   0,  0,  dphi3[1],   0,
                        0,  0,      dphi0[2],  0,  0,  dphi1[2],   0,  0,  dphi2[2],   0,  0,  dphi3[2],
                        dphi0[1],   dphi0[0],   0,  dphi1[1],   dphi1[0],   0, dphi2[1],   dphi2[0],   0, dphi3[1],   dphi3[0],   0,
                        0,  dphi0[2],   dphi0[1], 0,  dphi1[2],   dphi1[1], 0,  dphi2[2],   dphi2[1], 0,  dphi3[2],   dphi3[1],
                        dphi0[2],   0,  dphi0[0], dphi1[2],   0,  dphi1[0], dphi2[2],   0,  dphi2[0], dphi3[2],   0,  dphi3[0];
                
                
                
                //B*=(1.0/(6.0*m_V0));
                
                m_K.setZero();
                m_K = -m_V0*B.transpose()*C*B;
            }
            
            
            inline ~QuadratureExact() { }
            
            //integral rules for things that I want
            inline DataType getValue(const State<DataType> &state) {
                
                //Do nothing for now
                //returning the force which is really the negative gradient
                Eigen::Map<Eigen::VectorXd> v0 = mapDOFEigen(*m_qDofs[0], state);
                Eigen::Map<Eigen::VectorXd> v1 = mapDOFEigen(*m_qDofs[1], state);
                Eigen::Map<Eigen::VectorXd> v2 = mapDOFEigen(*m_qDofs[2], state);
                Eigen::Map<Eigen::VectorXd> v3 = mapDOFEigen(*m_qDofs[3], state);
                
                Eigen::Matrix<double, 12,1> q;
                q << v0, v1, v2, v3;

                
                return 0.5*q.transpose()*m_K*q;
                
            }
            
            template<typename Vector>
            inline void getGradient(Vector &f, const State<DataType> &state) {
                
                //returning the force which is really the negative gradient
                Eigen::Map<Eigen::VectorXd> v0 = mapDOFEigen(*m_qDofs[0], state);
                Eigen::Map<Eigen::VectorXd> v1 = mapDOFEigen(*m_qDofs[1], state);
                Eigen::Map<Eigen::VectorXd> v2 = mapDOFEigen(*m_qDofs[2], state);
                Eigen::Map<Eigen::VectorXd> v3 = mapDOFEigen(*m_qDofs[3], state);
                
                Eigen::Matrix<double, 12,1> q;
                q << v0, v1, v2, v3;
                Eigen::Matrix<double, 12,1> f0 = m_K*q;
                
                assign(f, f0, m_qDofs);
                
            }
            
            template<typename Matrix>
            inline void getHessian(Matrix &H, const State<DataType> &state) {
                
                assign(H, m_K, m_qDofs, m_qDofs);
                
            }
            
        protected:
            Eigen::Matrix<DataType, 12,12> m_K;
            double m_V0; //volume
            
        private:
        };
        
        //Body force due to gravity
        template<typename DataType, typename ShapeFunction>
        class QuadratureExact<DataType, BodyForceGravity<DataType, ShapeFunction> > :
        public BodyForceGravity<DataType, ShapeFunction> {
            
        public:
            
            using ShapeFunction::dphi;
            
            using BodyForceGravity<DataType, ShapeFunction>::m_rho;
            using BodyForceGravity<DataType, ShapeFunction>::m_g;
            using BodyForceGravity<DataType, ShapeFunction>::setBodyForceDensity;
            using ShapeFunction::m_qDofs;
            using ShapeFunction::m_qDotDofs;
            
            template<typename QDOFList, typename QDotDOFList>
            inline QuadratureExact(Eigen::MatrixXd &V, Eigen::MatrixXi &F, QDOFList &qDOFList, QDotDOFList &qDotDOFList) : BodyForceGravity<DataType, ShapeFunction>(V,F, qDOFList, qDotDOFList) {
                
                //need the volume of this tetrahedron
                Eigen::Matrix<DataType, 3,3> V0;
                V0 << Vert(1,0) - Vert(0,0), Vert(2,0) - Vert(0,0), Vert(3,0) - Vert(0,0),
                Vert(1,1) - Vert(0,1), Vert(2,1) - Vert(0,1), Vert(3,1) - Vert(0,1),
                Vert(1,2) - Vert(0,2), Vert(2,2) - Vert(0,2), Vert(3,2) - Vert(0,2);
                
                m_V0 = (1.0/6.0)*V0.determinant();
                
                assert(m_V0 > 0); //chewck tet is not inverted in reference config

                
            }
            
            
            inline ~QuadratureExact() { }
            
            //integral rules for things that I want
            template<typename DOFList>
            inline DataType getValue(const State<DataType> &state) {
                Eigen::Matrix<double, 12,1> f;
                getGradient(f, state);
                
                Eigen::Map<Eigen::VectorXd> v0 = mapDOFEigen(*m_qDofs[0], state);
                Eigen::Map<Eigen::VectorXd> v1 = mapDOFEigen(*m_qDofs[1], state);
                Eigen::Map<Eigen::VectorXd> v2 = mapDOFEigen(*m_qDofs[2], state);
                Eigen::Map<Eigen::VectorXd> v3 = mapDOFEigen(*m_qDofs[3], state);
                
                Eigen::Matrix<double, 12,1> q;
                q << v0, v1, v2, v3;
                
                return -q.transpose()*f;
            }
            
            template<typename Vector>
            inline void getGradient(Vector &f, const State<DataType> &state) {
                
                Eigen::Vector3d nodalF = m_rho*m_V0*0.25*m_g;
                Eigen::Matrix<double, 12,1> elementF = nodalF.replicate(4,1);
                //assign(f, elementF, std::array<DOFBase<DataType,0>, 4>{{*m_qDofs[0], *m_qDofs[1], *m_qDofs[2], *m_qDofs[3]}});
                assign(f, elementF, m_qDofs);
                
                
            }
            
            template<typename Matrix>
            inline void getHessian(Matrix &H, const State<DataType> &state) {
                
                
            }
            
        protected:
            double m_V0;
            
        private:
        };
        

    }
}

#endif /* QuadratureExact_h */
