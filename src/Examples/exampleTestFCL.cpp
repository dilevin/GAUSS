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
    
    //Read mesh to embed
    if(!igl::readOBJ(dataDir()+"/meshes/OBJ/garg.obj", V0, F0))
        std::cout<<"Failed to read embedded mesh \n";
    
    //Read mesh to embed
    if(!igl::readOBJ(dataDir()+"/meshes/OBJ/garg.obj", V1, F1))
        std::cout<<"Failed to read embedded mesh \n";
    
    //add to fcl BVH's using KDOPS
    fcl::BVHModel<fcl::KDOP<double, 24> > m0;
    fcl::BVHModel<fcl::KDOP<double, 24> > m1;
    
    //copy verts and
    m0.beginModel();
    m0.addSubModel(V0, F0);
    m0.endModel();
    
    m1.beginModel();
    m1.addSubModel(V1, F1);
    m1.endModel();
    
    fcl::Transform3<double> pose = fcl::Transform3<double>::Identity();
    
    //do collision detection
    fcl::CollisionRequest<double> request(1000, true);
    fcl::CollisionResult<double> result;
    int numContacts = fcl::collide(&m0,pose,&m1, pose, request, result);
    
    std::cout<<"Num Contacts: "<<numContacts<<"\n";
    
    std::vector<fcl::Contact<double>> contacts;
    result.getContacts(contacts);
    
    for(int i = 0; i < numContacts; ++i)
    {
        std::cout<<"Contact "<<i<<" "<<contacts[i].pos[0]<<" "<<contacts[i].pos[1]<<" "<<contacts[i].pos[2]<<"\n";
    }
    
    std::cout<<"Done \n";
    
}
