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
        
        template<typename DataType, typename Position1, typename Position2>
        class ForceSpringImpl
        {
        public:
            
            ForceSpringImpl(Position1 q0, Position2 q1, double l0, double k) : m_q0(0), m_q1(0) {
                m_q0 = q0;
                m_q1 = q1;
                m_l0 = l0;
                m_k = k;
            }
            
            //forces can give you the energy stored, the force itself and the hessian
            inline DataType getEnergy(State<DataType> &state) {
                
                //spring energy = k*(1.0 - l/l0).^2]
                double mag = (m_q1(state)-m_q0(state)).norm();
                //if(fabs(mag) < 1e-3) {
                  //  mag = 1e-3;
                //}
                
                DataType l = (1.0- mag/m_l0);
                return 0.5*m_k*l*l;
            }
            
            //forces always act on at least one DOF of the system this function returns which DOF the are acting on.
            inline auto & getDOF(unsigned int index) {
                
                if(index == 0) {
                    return *m_q0.getDOF();
                } else {
                    return *m_q1.getDOF();
                }
            }
            
            inline auto & getPosition0() {
                    return m_q0;
            }

            inline void setPosition0(Position1 q0) {
                    m_q0 = q0;
            }
            
            inline auto & getPosition1() {
                return m_q1;
            }

            inline void setPosition1(Position2 q1) {
                    m_q1 = q1;
            }
            
            inline unsigned int getNumDOF() { return 2; }
            
            const DataType & getStiffness() const { return m_k; }
            void setStiffness(const DataType &k) { m_k = k; }

            const DataType & getRestLength() const { return m_l0; }
            void setRestLength(const DataType &l0) { m_l0 = l0; }

            template <typename Vector>
            inline void getForce(Vector &f, State<DataType> &state) {
                
                auto q0 = m_q0(state);
                auto q1 = m_q1(state);
           
                double mag = (q1-q0).norm();
                
                
                if(fabs(mag) < 1e-8) {
                    mag = 1e-8;
               }
                
                
                Eigen::Vector3d fSpring = (m_k/m_l0)*(1.0/mag -1.0/m_l0)*(q1-q0);
                assign(f, fSpring, std::array<DOFBase<DataType,0> , 1>{{*m_q1.getDOF()}});
                fSpring = -fSpring;
                assign(f, fSpring, std::array<DOFBase<DataType,0> , 1>{{*m_q0.getDOF()}});
    
            }
            
            
            template <typename Matrix>
            inline void getStiffnessMatrix(Matrix &H, State<DataType> &state) {
                auto q0 = m_q0(state);
                auto q1 = m_q1(state);
                
                Eigen::Vector3x<DataType> dq = q1-q0;
               
                double mag = (q1-q0).norm();
                if(fabs(mag) < 1e-8) {
                    mag = 1e-8;
               }
                
                Eigen::Vector3x<DataType> dqN = dq/mag;
                
                //double strain = 1.0 - l/m_l0;
                double b = m_k/m_l0;
                double a = b*(1.0/mag - 1.0/m_l0);
                double c = b/(mag);
                
                Eigen::Matrix<double, 6,3> B;
                B << -1,  0,  0,
                      0, -1,  0,
                      0,  0, -1,
                      1,  0,  0,
                      0,  1,  0,
                      0,  0,  1;
                
                Eigen::Matrix<double, 6,6> Hspring;
                Hspring = a*B*B.transpose() - c*B*dqN*dqN.transpose()*B.transpose();
                assign(H, Hspring, std::array<DOFBase<DataType,0> , 1>{{*m_q0.getDOF()}}, std::array<DOFBase<DataType,0> , 1>{{*m_q1.getDOF()}});
                
            }
            
        protected:
            
            DataType m_l0, m_k; //spring original length and stiffness
            Position1 m_q0; //spring end points
            Position2 m_q1;
            
        private:
        };
        
        template<typename DataType, typename Position1, typename Position2>
        using ForceSpring = Force<DataType, ForceSpringImpl<DataType, Position1, Position2> >;
        
        template<typename DataType>
        using PosParticle = PositionEigen<DataType, DOFParticle<DataType> >;
        
        template<typename DataType>
        using ForceSpringParticles = Force<DataType, ForceSpringImpl<DataType,
                                                           PosParticle<DataType>,
                                                           PosParticle<DataType> > >;
    }
}

#endif /* ForceSpring_h */
