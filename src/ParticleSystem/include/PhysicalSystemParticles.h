#ifndef PHYSICALSYSTEMPARTICLES_H
#define PHYSICALSYSTEMPARTICLES_H

#include "PhysicalSystem.h"
#include "DOFParticle.h"
#include <array>

//Simple test which is a single particle (just to shake out bugs in the design
namespace Gauss {
    namespace ParticleSystem {
        
        template<typename DataType>
        class PhysicalSystemParticleSingleImpl
        {
        public:
            
            //temporary global indices until I update the state to give these to me
            //automatically
            PhysicalSystemParticleSingleImpl() : m_x(0), m_xDot(0) {
                m_mass = 1.5; //temporary set mass to avoid weird failures
            }
            
            template<typename Assembler>
            inline void getMassMatrix(Assembler &assembler, const State<DataType> &state) const {
                assign(assembler, m_mass, std::array<DOFBase<DataType,1>, 1>{{m_xDot}}, std::array<DOFBase<DataType,1>, 1>{{m_xDot}});
            }
            
            template<typename Assembler>
            inline void getStiffnessMatrix(Assembler &assembler, const State<DataType> &state) const {
                
                //DO NOTHING ... a single particle has no stiffness
            }
            
            template<typename Assembler>
            inline void getForce(Assembler &assembler, const State<DataType> &state) const {
                
                //DO NOTHING ... a single particle has no stiffness
            }
            
            const double & getMass() const { return m_mass; }
            void setMass(const double &mass) { m_mass = mass; }
            
            DOFParticle<DataType,0> & getQ() { return m_x; }
            const DOFParticle<DataType,0> & getQ() const { return m_x; }
            
            DOFParticle<DataType,1> & getQDot() { return m_xDot; }
            const DOFParticle<DataType,1> & getQDot() const { return m_xDot; }
            
        protected:
            
            DataType m_mass; //mass of particle
            DOFParticle<DataType,0> m_x;
            DOFParticle<DataType,1> m_xDot;
            
        private:
        };
        
        template<typename DataType>
        using  PhysicalSystemParticleSingle = PhysicalSystem<DataType, PhysicalSystemParticleSingleImpl<DataType> >;

    }
}


#endif
