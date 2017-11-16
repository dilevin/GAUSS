//
//  ForceSpring.h
//  Gauss
//
//  Created by David Levin on 3/4/17.
//
//

#ifndef ForceGravity_h
#define ForceSpring_h

#include <ForceExternal.h>
#include <DOFParticle.h>
#include <UtilitiesEigen.h>
#include <unsupported/Eigen/AutoDiff>

namespace Gauss {
    namespace ParticleSystem {
        
        template<typename DataType>
        class ForceGravityImpl
        {
        public:
            
            ForceGravityImpl(ParticleSystem::DOFParticle<DataType> *q0, DataType mass, Eigen::Vector3x<DataType> g)  {
                m_q0 = q0;
                m_g = mass*g;
            }
            
            //forces can give you the energy stored, the force itself and the hessian
            template<typename Scalar>
            inline void getEnergy(Scalar &e,  State<DataType> &state) {
                
            }
            
            //forces always act on at least one DOF of the system this function returns which DOF the are acting on.
            inline auto & getDOF(unsigned int index) {
                
                return *m_q0;
            }
            
            inline unsigned int getNumDOF() { return 1; }
            
            template <typename Vector>
            inline void getForce(Vector &f, State<DataType> &state) {
                assign(f, m_g, std::array<DOFBase<DataType,0> , 1>{{*m_q0}});
            }
            
            
            template <typename Matrix>
            inline void getStiffnessMatrix(Matrix &H, State<DataType> &state) {
                
                
            }
            
        protected:
            
            ParticleSystem::DOFParticle<DataType> *m_q0;
            Eigen::Vector3x<DataType> m_g;
            
        private:
        };
        
        
        template<typename DataType>
        using ForceParticlesGravity = Force<DataType, ForceGravityImpl<DataType> >;
    }
}

#endif /* ForceSpring_h */
