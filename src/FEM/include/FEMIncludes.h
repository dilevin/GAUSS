//
//  FEMIncludes.h
//  Gauss
//
//  Created by David Levin on 4/8/17.
//
//

#ifndef FEMIncludes_h
#define FEMIncludes_h

#include <DOFParticle.h>
#include <PhysicalSystemFEM.h>
#include <Element.h>
#include <Energy.h>
#include <EnergyNeohookean.h>

//Quadrature rules
#include <Quadrature.h>
#include <QuadratureExact.h>
#include <QuadratureTetConstant.h>
#include <QuadratureHex8.h>
#include <QuadraturePlaneTri1.h>

#include <ShapeFunction.h>
#include <ShapeFunctionHexTrilinear.h>
#include <ShapeFunctionPlaneLinear.h>
#include <ElementTypes.h>
#include <UtilitiesIO.h> // I want to be able to load tetgen files

//Utils
#include <UtilitiesFEM.h>

#endif /* FEMIncludes_h */
