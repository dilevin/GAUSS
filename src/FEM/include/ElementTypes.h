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
        using NeohookeanTet = ElementBase<DataType, 4, QuadratureExact, QuadratureTetConstant, QuadratureTetConstant, EnergyKineticNonLumped, EnergyNeohookean, BodyForceGravity, ShapeFunctionLinearTet>;

        template<typename DataType>

        using NeohookeanHFixedTet = ElementBase<DataType, 4, QuadratureExact, QuadratureTetConstant, QuadratureTetConstant, EnergyKineticNonLumped, EnergyNeohookeanHFixed, BodyForceGravity, ShapeFunctionLinearTet>;

        template<typename DataType>
        using MuscleTet = ElementBase<DataType, 4, QuadratureExact, QuadratureTetConstant, QuadratureTetConstant, EnergyKineticNonLumped, EnergyMuscle, BodyForceGravity, ShapeFunctionLinearTet>;

        template<typename DataType>
        using StvkTet = ElementBase<DataType, 4, QuadratureExact, QuadratureTetConstant, QuadratureTetConstant, EnergyKineticNonLumped, EnergyStvk, BodyForceGravity, ShapeFunctionLinearTet>;
        
        template<typename DataType>
        using LinearHex = Element<DataType, 8, QuadratureHex8, EnergyKineticNonLumped, EnergyLinearElasticity, BodyForceGravity, ShapeFunctionHexTrilinear>;

        template<typename DataType>
        using NeohookeanHex = Element<DataType, 8, QuadratureHex8, EnergyKineticNonLumped, EnergyNeohookean, BodyForceGravity, ShapeFunctionHexTrilinear>;
        
        template<typename DataType>
        using LinearPlaneStrainTri = ElementBase<DataType, 3, QuadratureExact, QuadraturePlaneTri1, QuadraturePlaneTri1, EnergyKineticNonLumped, EnergyLinearElasticity, BodyForceGravity, ShapeFunctionPlaneLinear>;

        template<typename DataType>
        using MuscleTri = ElementBase<DataType, 3, QuadratureExact, QuadraturePlaneTri1, QuadraturePlaneTri1, EnergyKineticNonLumped, EnergyMuscle, BodyForceGravity, ShapeFunctionPlaneLinear>;
        
        template<typename DataType, template<typename A, typename B> class EnergyPS>
        using FEMPrincipalStretchTet = ElementBase<DataType, 4, QuadratureExact, QuadratureTetConstant, QuadratureTetConstant, EnergyKineticNonLumped, EnergyPS, BodyForceGravity, ShapeFunctionLinearTet>;

    }
}
#endif /* ElementTypes_h */
