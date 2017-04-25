//
//  ForceSpring.h
//  Gauss
//
//  Created by David Levin on 3/4/17.
//
//

#ifndef ForceSpring_h
#define ForceSpring_h

#include <ForceExternal.h>
#include <DOFParticle.h>
#include <UtilitiesEigen.h>
#include <unsupported/Eigen/AutoDiff>

namespace Gauss {
    namespace ParticleSystem {
        
        template<typename DataType> 
        class ForceSpringImpl
        {
        public:
            
            ForceSpringImpl(DOFParticle<DataType> *q0, DOFParticle<DataType> *q1, double l0, double k) {
                m_q[0] = q0;
                m_q[1] = q1;
                m_l0 = l0;
                m_k = k;
            }
            
            //forces can give you the energy stored, the force itself and the hessian
            template<typename Scalar>
            inline void getEnergy(Scalar &e,  State<DataType> &state) {
                Eigen::Map<Eigen::VectorXd> q0 = mapDOFEigen(*m_q[0], state);
                Eigen::Map<Eigen::VectorXd> q1 = mapDOFEigen(*m_q[1], state);
                
                //spring energy = k*(1.0 - l/l0).^2]
                DataType l = (1.0-(q1-q0).norm());
                e = 0.5*m_k*l*l;
            }
            
            //forces always act on at least one DOF of the system this function returns which DOF the are acting on.
            inline auto & getDOF(unsigned int index) {
                
                return *(m_q[index]);
            }
            
            inline unsigned int getNumDOF() { return 2; }
            
            template <typename Vector>
            inline void getForce(Vector &f, State<DataType> &state) {
                
                Eigen::Map<Eigen::VectorXd> q0 = mapDOFEigen(*m_q[0], state);
                Eigen::Map<Eigen::VectorXd> q1 = mapDOFEigen(*m_q[1], state);
                
                double l = (q1-q0).norm();
                double strain = 1.0 - l;
                Eigen::Vector3d fSpring = (m_k/m_l0)*strain*(q1-q0)/l;
                assign(f, fSpring, std::array<DOFBase<DataType,0> , 1>{{*m_q[1]}});
                fSpring = -fSpring;
                assign(f, fSpring, std::array<DOFBase<DataType,0> , 1>{{*m_q[0]}});
    
            }
            
            
            template <typename Matrix>
            inline void getStiffnessMatrix(Matrix &H, State<DataType> &state) {
                
            }
            

            
        protected:
            
            DataType m_l0, m_k; //spring original length and stiffness
            DOFParticle<DataType> *m_q[2]; //spring end points
            
        private:
        };
        
        template<typename DataType>
        using ForceSpring = Force<DataType, ForceSpringImpl<DataType> >;
    }
}

#endif /* ForceSpring_h */
