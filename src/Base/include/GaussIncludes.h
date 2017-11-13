//
//  GaussIncludes.h
//  Gauss
//
//  Created by David Levin on 4/8/17.
//
//

#ifndef GaussIncludes_h
#define GaussIncludes_h

#include <iostream>
#include <PhysicalSystem.h>
#include <MultiVector.h>
#include <World.h>
#include <Assembler.h>
#include <Utilities.h>
#include <UtilitiesEigen.h>
#include <UtilitiesGeometry.h>
#include <UtilitiesBase.h>
#include <UtilitiesMATLAB.h>

#ifdef GAUSS_PARDISO
#include <SolverPardiso.h>
#endif

#include <AssemblerParallel.h>


#endif /* GaussIncludes_h */
