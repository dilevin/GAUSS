#ifndef COLLISIONSFCL_H
#define COLLISIONSFCL_H

//FCL
#include "fcl/math/bv/utility.h"
#include "fcl/narrowphase/collision.h"
#include "fcl/narrowphase/detail/gjk_solver_indep.h"
#include "fcl/narrowphase/detail/gjk_solver_libccd.h"
#include "fcl/narrowphase/detail/traversal/collision_node.h"

//GAUSS Stuff
#include <GaussIncludes.h>
#include <UtilitiesContact.h>
#define MAX_CONTACTS 1000

// Collision detector fo    r triangle meshes using the flexibile collision library (https://github.com/flexible-collision-library/fcl)
namespace Gauss {
    namespace Collisions {
        
        //Some useful methods for plane-object collision checking
        template<typename DataType>
        class CollisionFCLImpl {
            
        public:
            
            template<typename World>
            CollisionFCLImpl(World &world) {
                
                //scan through the world and put all the object into FCL collision structures
                //std::vector< fcl::BVHModel<fcl::OBBRSS<DataType> > *> &bvhList = m_bvhList;
                //std::vector< SystemIndex> &indexList = m_indexList;
                m_firstTime = true;
                
                             
            }
            
            ~CollisionFCLImpl() {
                
                //clean up FCL data structures
            }
            
            template<typename World>
            void detectCollisions(World &world);
            
            inline unsigned int getNumCollisions() const { return m_sharedList.size(); }
            
            inline auto & getCollisionsObjectA() { return m_objAList; }
            inline auto & getCollisionsObjectB() { return m_objBList; }
            inline auto & getSharedInfo() { return m_sharedList; }
            
            inline const auto & getCollisionsObjectA() const { return m_objAList; }
            inline const auto & getCollisionsObjectB() const { return m_objBList; }
            inline const auto & getSharedInfo() const { return m_sharedList; }
            
            
        protected:
            
            //only vertex collisions for now
            MultiVector<ObjectCollisionInfo<DataType, 1> > m_objAList; //collision info for object A
            MultiVector<ObjectCollisionInfo<DataType, 1> > m_objBList; //collision info for object B
            std::vector<SharedCollisionInfo<DataType>> m_sharedList; //normals and world space positions;
            
            //FCL collision data structures
            std::vector< fcl::BVHModel<fcl::OBBRSS<DataType> > > m_bvhList;
            std::vector< SystemIndex > m_indexList; 
            
            bool m_firstTime;
            
        private:
            
        };
    }
}

template<typename DataType>
template<typename World>
void Gauss::Collisions::CollisionFCLImpl<DataType>::detectCollisions(World &world) {
    
    m_sharedList.clear();
    m_objAList.clear();
    m_objBList.clear();
    
    //MultiVector<Gauss::Collisions::ObjectCollisionInfo<DataType, 0> > &objAList = m_objAList;
    //MultiVector<Gauss::Collisions::ObjectCollisionInfo<DataType, 0> > &objBList = m_objBList;
    //std::vector<Gauss::Collisions::SharedCollisionInfo<DataType> > &sharedList = m_sharedList;
    
    std::vector< fcl::BVHModel<fcl::OBBRSS<DataType> > > &bvhList = m_bvhList;
    std::vector< SystemIndex > &indexList = m_indexList;
    
    if(m_firstTime) {
        
        m_bvhList.clear();
        m_indexList.clear();
        
        forEachIndex(world.getSystemList(), [&world, &bvhList, &indexList](auto type, auto index, auto &a) {
            
            Eigen::MatrixXx<DataType> positions;
            positions.resize(a->getGeometry().first.rows(),3);
            
            for(unsigned int ii=0; ii<a->getGeometry().first.rows(); ++ii) {
                positions.row(ii) = a->getPosition(world.getState(),ii);
            }
            
            bvhList.push_back(fcl::BVHModel<fcl::OBBRSS<DataType> >());
            
            bvhList[bvhList.size()-1].beginModel();
            bvhList[bvhList.size()-1].addSubModel(positions, a->getGeometry().second);
            bvhList[bvhList.size()-1].endModel();
            
            
            //bvhList.push_back(bvh);
            indexList.push_back(SystemIndex(type, index));
        });
        
        m_firstTime = false;
    }
    
    //update bvhs
     unsigned int bvhIndex = 0;
    //fill in collision data using retrived indices and collision data
    forEachIndex(world.getSystemList(), [&bvhList, &bvhIndex, &world](auto type, auto index, auto &a) {
        
        bvhList[bvhIndex].beginUpdateModel();
        for(unsigned int ii=0; ii<a->getGeometry().first.rows(); ++ii) {
            bvhList[bvhIndex].updateVertex(a->getPosition(world.getState(), ii));
        }
        
        bvhList[bvhIndex].endUpdateModel(true, true);
        
        bvhIndex++;
            
    });
    
    //n^2 collision detection
    fcl::Transform3<DataType> pose0 = fcl::Transform3<DataType>::Identity();
    fcl::Transform3<DataType> pose1 = fcl::Transform3<DataType>::Identity();
    fcl::CollisionRequest<DataType> request(MAX_CONTACTS,true);
    fcl::CollisionResult<DataType> result;
    
    for(unsigned int obj0=0; obj0<m_bvhList.size(); ++obj0) {
        for(unsigned int obj1=obj0+1; obj1<m_bvhList.size(); ++obj1) {
            
            //do collision detection
            fcl::collide(&m_bvhList[obj0],pose0, &m_bvhList[obj1], pose1, request, result);
            
            int numContacts = result.numContacts();
            
            for(unsigned int ii=0; ii<numContacts; ++ii) {
                
            //add each new contact to the contact list here
            //b1, b2 = collision primitive in obj1 and obj2
            //normal = contact normal
            //x = world space contact position
                
                //std::cout<<"Position: \n"<<result.getContact(ii).pos<<", \nNormal: \n"<<result.getContact(ii).normal<<"\n";
                m_sharedList.push_back(SharedCollisionInfo<DataType>(-result.getContact(ii).normal,result.getContact(ii).pos));
                m_objAList.add(ObjectCollisionInfo<DataType,1>(m_sharedList.size()-1,m_indexList[obj0], result.getContact(ii).b1));
                m_objBList.add(ObjectCollisionInfo<DataType,1>(m_sharedList.size()-1,m_indexList[obj1], result.getContact(ii).b2));
            }
        }
    }
    
    
    
}

#endif
