#ifndef DOFPARTICLE_H
#define DOFPARTICLE_H

#include "DOF.h"
#include "Utilities.h"

namespace Gauss {
    namespace ParticleSystem {
        
        /**
         Implementation for particle in R3 (can also be used for particle velocity)
         */
        template<typename DataType, unsigned int Property=0>
        class DOFParticleImpl
        {
        public:
    
            explicit DOFParticleImpl() {
                
            }
        
            
            //The t parameter is optional and is used for spacetime degrees of freedom
            inline unsigned int getNumScalarDOF() const {
                return 3;
            }
            
            //get dataPtr for this DOF
            //inline std::tuple<DataType *, unsigned int> getPtr(State<DataType> &state, unsigned int offset) {
             //   return state. template getStatePtr<Property>(offset);
            //}
            
        protected:
        

        private:
        
        };
        
        template<typename DataType, unsigned int PropertyIndex=0>
        using DOFParticle = DOF<DataType, DOFParticleImpl, PropertyIndex>;
    }
}
#endif
