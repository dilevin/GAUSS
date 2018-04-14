//
//  ConstraintContact.h
//  Gauss
//
//  Created by David Levin on 11/18/17.
//
//

#ifndef ConstraintContact_h
#define ConstraintContact_h

#include <UtilitiesEigen.h>

namespace Gauss {
    namespace Collisions {
   
        template<typename DataType>
        struct SharedCollisionInfo {
            
            inline SharedCollisionInfo() {
                
            }
            
            inline SharedCollisionInfo(const Eigen::Vector3x<DataType> &normal, const Eigen::Vector3x<DataType> position) {
                m_normal = normal;
                m_contactPoint = position;
            }
            
            inline const Eigen::Vector3x<DataType> & getNormal() const { return m_normal; }
            inline const Eigen::Vector3x<DataType> & getPosition() const { return m_contactPoint; }
            
            //important data (I borrowed the idea of this stucture from the flexible-collision-library)
            Eigen::Vector3x<DataType> m_normal; //as a rule m_normal is pointing away from object A
            Eigen::Vector3x<DataType> m_contactPoint;

        };
        
        //Store a collision pair, gets filled in by collision detector
        //Type is the type of collision for this object
        //0 - vertex collision
        //1 - element collision
        //2 - position only collision
        template<typename DataType, unsigned int Type>
        struct ObjectCollisionInfo {
            
            static constexpr unsigned int collisionType = Type;
            
            inline ObjectCollisionInfo() {
                m_indices[0] = -1;
                m_indices[1] = -1;
                m_infoIndex = -1;
            }
            
            inline ObjectCollisionInfo(int infoIndex, SystemIndex A, int geometryIndex) : m_infoIndex(infoIndex) {
                m_A = A;
                m_indices[0] = -1;
                m_indices[1] = -1;
                m_indices[Type] = geometryIndex;
                m_infoIndex = infoIndex;

            }
            
            inline const SystemIndex & getObject() const { return m_A; }
            inline const int & getData(unsigned int i) const { assert(i < 2); return m_indices[i]; }
            inline const int & getShared() const { return m_infoIndex; }
            SystemIndex m_A;
            
            //additional information
            //0-contact element object
            //1- vertex id
            
            //These are empty of object is not a mesh geometry
            std::array<int, 3> m_indices; //extra int for overflow to avoid an if in the constructor
            int m_infoIndex;
            
        };
        
    }
}

#endif /* ConstraintContact_h */
