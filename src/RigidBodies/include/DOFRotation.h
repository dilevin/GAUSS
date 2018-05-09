//
//  DOFRotation.h
//  Gauss
//
//  Created by David Levin on 5/8/18.
//

#ifndef DOFRotation_h
#define DOFRotation_h

#include "DOF.h"
#include "Utilities.h"

namespace Gauss {
    namespace RigidBodies {
        
        /**
         Implementation Rotation degree of freedom stored as a four value quaternion
         */
        template<typename DataType, unsigned int Property=0>
        class DOFRotationImpl
        {
        public:
            
            explicit DOFRotationImpl() {
                
            }
            
            
            //The t parameter is optional and is used for spacetime degrees of freedom
            inline unsigned int getNumScalarDOF() const {
                return 4;
            }
            
            //get dataPtr for this DOF
            //inline std::tuple<DataType *, unsigned int> getPtr(State<DataType> &state, unsigned int offset) {
            //   return state. template getStatePtr<Property>(offset);
            //}
            
        protected:
            
            
        private:
            
        };
        
        template<typename DataType, unsigned int PropertyIndex=0>
        using DOFRotation = DOF<DataType, DOFRotationImpl, PropertyIndex>;
    }
}

#endif /* DOFRotation_h */
