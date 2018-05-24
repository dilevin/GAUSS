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
    
    /**
     Implementation Rotation degree of freedom stored as a four value quaternion
     */
    template<typename DataType, unsigned int Property=0>
    class DOFRotationImpl
    {
    public:
        
        explicit DOFRotationImpl() {
            
        }
        
        inline unsigned int getNumScalarDOF() const {
            return 4;
        }
        
    protected:
        
        
    private:
        
    };
    
    template<typename DataType, unsigned int PropertyIndex=0>
    using DOFRotation = DOF<DataType, DOFRotationImpl, PropertyIndex>;
}

#endif /* DOFRotation_h */
