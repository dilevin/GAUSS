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

    }
}
#endif /* ElementTypes_h */
