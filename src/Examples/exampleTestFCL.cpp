#include <iostream>

#include <functional>

#include <Qt3DIncludes.h>
#include <GaussIncludes.h>
#include <FEMIncludes.h>

//Any extra things I need such as constraints
#include <ConstraintFixedPoint.h>
#include <TimeStepperEulerImplicitLinear.h>
#include <Embedding.h>

using namespace Gauss;
using namespace FEM;
using namespace ParticleSystem; //For Force Spring


//FCL
#include "fcl/math/bv/utility.h"
#include "fcl/narrowphase/collision.h"
#include "fcl/narrowphase/detail/gjk_solver_indep.h"
#include "fcl/narrowphase/detail/gjk_solver_libccd.h"
#include "fcl/narrowphase/detail/traversal/collision_node.h"

int main(int argc, char **argv) {

    //load two tet meshes
    Eigen::MatrixXd V0, V1;
    Eigen::MatrixXi F0, F1;
    
    readTetgen(V0, F0, dataDir()+"/meshesTetgen/Gargoyle/gargNested2.node", dataDir()+"/meshesTetgen/Gargoyle/gargNested2.ele");
    
    readTetgen(V1, F1, dataDir()+"/meshesTetgen/Gargoyle/gargNested2.node", dataDir()+"/meshesTetgen/Gargoyle/gargNested2.ele");
    
    //add to fcl BVH's using KDOPS
    
    //do collision detection
    
    
}
