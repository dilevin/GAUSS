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
        using NeohookeanTet = ElementBase<DataType, 4, QuadratureExact, QuadratureTetConstant, EnergyKineticNonLumped, EnergyNeohookean, BodyForceGravity, ShapeFunctionLinearTet>;

        
        template<typename DataType>
        using LinearHex = Element<DataType, 8, QuadratureHex8, EnergyKineticNonLumped, EnergyLinearElasticity, BodyForceGravity, ShapeFunctionHexTrilinear>;

        template<typename DataType>
        using NeohookeanHex = Element<DataType, 8, QuadratureHex8, EnergyKineticNonLumped, EnergyNeohookean, BodyForceGravity, ShapeFunctionHexTrilinear>;
        
        template<typename DataType>
        using LinearPlaneStrainTri = ElementBase<DataType, 3, QuadratureExact, QuadraturePlaneTri1, EnergyKineticNonLumped, EnergyLinearElasticity, BodyForceGravity, ShapeFunctionPlaneLinear>;

    }
}
#endif /* ElementTypes_h */
