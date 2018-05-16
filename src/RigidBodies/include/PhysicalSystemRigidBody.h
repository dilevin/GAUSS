//
//  PhysicalSystemRigidBody.h
//  Gauss
//
//  Created by David Levin on 5/8/18.
//

#ifndef PhysicalSystemRigidBody_h
#define PhysicalSystemRigidBody_h

//A dynamic rigid body
#include <vector>
#include <DOFParticle.h>
#include <DOFRotation.h>
#include <DOFPair.h>
#include <UtilitiesEigen.h>
#include <UtilitiesRigidBodies.h>

namespace Gauss {
    namespace RigidBodies {
    
        template<typename DataType>
        class PhysicalSystemRigidBodyImpl
        {
            
        public:
            
            using GammaMatrix = Eigen::Matrix<DataType, 3,6>;
            
            //temporary global indices until I update the state to give these to me
            //automatically
            PhysicalSystemRigidBodyImpl(const Eigen::Ref<Eigen::MatrixXd> &V, const Eigen::Ref<Eigen::MatrixXi> &F, DataType density=1.0) {
                
                m_V = V;
                m_F = F;
                m_numVerts = m_V.rows();
                m_numFaces = m_F.rows();
                
                assert(m_V.cols() == 3); //3D only for now
                
                Eigen::Vector3x<DataType> inertia;
                m_mass = computeMoments(m_com, inertia, m_R0, m_V, m_F);
                m_mass *= density;
                inertia *= density;
                m_massMatrix.setConstant(m_mass);
                m_massMatrix.segment(0,3) = inertia;
                
                //subtract center of mass off of vertex positions
                #pragma omp parallel for
                for(unsigned int ii = 0; ii < V.rows(); ++ii)
                {
                    m_V.row(ii) -= m_com.transpose();
                }
                
                //fine I'll rotate the stupid mesh
                m_V = (m_R0.transpose()*m_V.transpose()).transpose();
                
            }
            
            ~PhysicalSystemRigidBodyImpl() {
                
            }
            
            DataType getEnergy(const State<DataType> &state) const {
                
            }
            
            DataType getStrainEnergy(const State<DataType> &state) const {
                
                return 0.0;
            }
            
            template<typename Assembler>
            inline void getMassMatrix(Assembler &assembler, const State<DataType> &state) const {
                
                //Rigid body mass matrix
                assign(assembler, m_massMatrix.asDiagonal().toDenseMatrix().eval(), getQDot(0), getQDot(0));
            }
            
            //the rigid body jacobian
            
            template<typename Assembler>
            inline void getStiffnessMatrix(Assembler &assembler, const State<DataType> &state) const {
                
                //do nothing, rigid body has no constitutive model
            }
            
            template<typename Assembler>
            inline void getForce(Assembler &assembler, const State<DataType> &state) const {
                
                getBodyForce(assembler, state);
            }
            
            template<typename Assembler>
            inline void getInternalForce(Assembler &assembler, const State<DataType> &state) const {
                
                //centripedal force and coriolis force
            }
            
            template<typename Assembler>
            inline void getBodyForce(Assembler &assembler, const State<DataType> &state) const {
                
                //gravity goes here
                Eigen::Matrix<DataType, 6,1> g;
                g << 0,0,0, 0, m_mass*(-9.8), 0.0;
                g.segment(3,3) = mapDOFEigenQuat(m_q.first(), state).toRotationMatrix().transpose()*m_R0.transpose()*g.segment(3,3);
                std::cout<<"W: "<<mapDOFEigenQuat(m_q.first(), state).x()<<" "<<mapDOFEigenQuat(m_q.first(), state).y()<<" "<<mapDOFEigenQuat(m_q.first(), state).z()<<" "<<mapDOFEigenQuat(m_q.first(), state).w()<<"\n";
                std::cout<<"ROTATION MATRIX: \n"<<mapDOFEigenQuat(m_q.first(), state).toRotationMatrix()<<"\n";
                assign(assembler, g, getQDot(0));
            }
            
            //Degree-of-freedom access
            inline auto & getQ() { return m_q; }
            inline const auto & getQ() const { return m_q; }
            
            inline auto & getQDot() { return m_qDot; }
            inline const auto & getQDot() const { return m_qDot; }
            
            //get function supporting a vertex (these return arrays in order to slot directly into assemblers)
            inline decltype(auto) getQ(unsigned int vertexId) const {
                std::array<const DOFBase<DataType,0> *,1> toReturn = {{&m_q}};
                return toReturn;
            }
            
            inline decltype(auto) getQDot(unsigned int vertexId) const {
                std::array<const DOFBase<DataType,1> *,1> toReturn = {{&m_qDot}};
                return toReturn;
            }
            
            
            template<typename Vector>
            inline decltype(auto) getQ(Vector &x, unsigned int elementId) const {
                std::cout<<"Error not implemented \n";
                exit(0);
                std::array<const DOFBase<DataType,0> *, 1> toReturn = {{&m_q[elementId]}};
                return toReturn;
            }
            
            template<typename Vector>
            inline decltype(auto) getQDot(Vector &x, unsigned int elementId) const {
                std::cout<<"Error not implemented \n";
                exit(0);
                std::array<const DOFBase<DataType,1> *,1> toReturn = {{&m_qDot[elementId]}};
                return toReturn;
            }
            
            //Geometry
            inline const auto getPosition(const State<DataType> &state, unsigned int vertexId) const {
                return  mapDOFEigen(m_q.first(), state).toRotationMatrix()*m_R0*(m_V.row(vertexId).transpose() + mapDOFEigen(m_q.second(), state));
            }
            
            inline const auto getVelocity(const State<DataType> &state, unsigned int vertexId) const {
                return getDPDQ(state,vertexId)*mapDOFEigen(m_qDot, state);
            }
            
            inline const auto getDPDQ(const State<DataType> &state, unsigned int vertexId) const {
                
                GammaMatrix gamma;
                gamma.setZero();
                
                gamma.block(0,0,3,3) << 0, m_V(vertexId,2), -m_V(vertexId,1), -m_V(vertexId,2), 0, m_V(vertexId,0), m_V(vertexId,1), -m_V(vertexId,0), 0;
                gamma.block(0,3,3,3).setIdentity();
                return mapDOFEigen(m_q.first(), state).toRotationMatrix()*m_R0*gamma;
            }
            
            //Rigid body Jacobian goes here
            inline const auto getDPDQ(const State<DataType> &state, unsigned int elementId, const Eigen::Vector3x<DataType> &pos) const {
                
                std::cout<<"position-based DPDQ not implemented in rigid body system \n";
                exit(0);
            }
            
            inline auto getGeometry() { return std::make_pair(std::ref(m_V), std::ref(m_F)); }
            
        protected:
            
            //Mesh
            unsigned int m_numVerts, m_numFaces;
            Eigen::MatrixXd m_V;
            Eigen::MatrixXi m_F;
            
            DataType m_mass;
            Eigen::Matrix33x<DataType> m_R0; //initial rotation from inertia frame to initial state
            Eigen::Vector3x<DataType> m_com; //center of mass in body frame
            Eigen::Vector6x<DataType> m_massMatrix; //mass matrix is diagonal for rigid bodies (yay!)
            
            //positions are a rotation and a particle (for the translation) -- stored in body frame
            DOFPair<DataType, DOFRotation, DOFParticle, 0> m_q;
            
            //velocities are an angular velocity (particle) and a linear velocity (particle) -- stored in body frame
            DOFPair<DataType, DOFParticle, DOFParticle, 1> m_qDot;
            
            
        private:
            
        };
        
        template<typename DataType>
        using  PhysicalSystemRigidBody = PhysicalSystem<DataType, PhysicalSystemRigidBodyImpl<DataType> >;
        
    }
}

#endif /* PhysicalSystemRigidBody_h */
