//
//  ElementTypes.h
//  Gauss
//
//  Created by David Levin on 4/8/17.
//
//

#ifndef ElementTypes_h
#define ElementTypes_h

namespace Gauss {
    namespace FEM {
        template<typename DataType>
        using LinearTet = Element<DataType, 4, QuadratureExact, EnergyKineticNonLumped, EnergyLinearElasticity, BodyForceGravity, ShapeFunctionLinearTet>;
        
        template<typename DataType>
        using LinearHex = Element<DataType, 8, QuadratureNone, EnergyKineticNonLumped, EnergyLinearElasticity, BodyForceGravity, ShapeFunctionHexTrilinear>;

    }
}
#endif /* ElementTypes_h */
