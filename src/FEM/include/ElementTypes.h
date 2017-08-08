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
        using LinearHex = Element<DataType, 8, QuadratureHex8, EnergyKineticNonLumped, EnergyLinearElasticity, BodyForceGravity, ShapeFunctionHexTrilinear>;

        template<typename DataType>
        using NeohookeanHex = Element<DataType, 8, QuadratureHex8, EnergyKineticNonLumped, EnergyNeohookean, BodyForceGravity, ShapeFunctionHexTrilinear>;

    }
}
#endif /* ElementTypes_h */
