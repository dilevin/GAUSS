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
   
        struct SystemIndex {
            
            inline SystemIndex() {
                m_type = -1;
                m_index = 0;
            }
            
            inline SystemIndex(unsigned int type, unsigned int index) {
                m_type = type;
                m_index = index;
            }
            
            inline int getType() const { return m_type; }
            inline unsigned int getIndex() const { return m_index; }
            
            inline int & getType() { return m_type; }
            inline unsigned int & getIndex() { return m_index; }
            
            int m_type; //-1 is fixed object, doesn't need collision response
            int m_index; //index of object in respective systems list
            
            
        };

        struct CollisionInfo {

            inline CollisionInfo(SystemIndex system, unsigned int elementId = 0, unsigned int vertexId = 0);
            
            inline unsigned int
            SystemIndex m_system;
            unsigned int m_vertexId;
            unsigned int m_elementId;
            
        };

        //Store a collision pair, gets filled in by collision detector 
        template<typename DataType>
        struct CollisionPair {

            inline CollisionPair() {
                
            }
            
            inline CollisionPair(SystemIndex A, SystemIndex B, Eigen::Vector3x normal, Eigen::Vector3x point,
                                  int elementA = -1,  int elementB = -1, int vertexA = -1, int vertexB = -1) {
                m_A = A;
                m_B = B;
                m_normal = normal;
                m_contactPoint = contactPoint
                m_indices[0] = elementA;
                m_indices[1] = elementB;
                m_indices[2] = vertexA;
                m_indices[3] = vertexB;
            }
            
            inline const SystemIndex & getA() const { return m_A; }
            inline const SystemIndex & getB() const { return m_B; }
            inline const Eigen::Vector3x<DataType> & getNormal() const { return m_normal; }
            inline const Eigen::Vector3x<DataType> & getPoint() const { return m_constactPoint; }
            
            SystemIndex m_A;
            System/index m_B;
            
            //important data (I borrowed the idea of this stucture from the flexible-collision-library)
            Eigen::Vector3x<DataType> m_normal;
            Eigen::Vector3x<DataType> m_contactPoint;
            
            //additional information
            //0-contact element object A (if has mesh)
            //1- element for B
            //2-contact vertex object A (if has mesh)
            //3- vertex for B
            
            //These are empty of object is not a mesh geometry
            std::array<int, 4> m_indices;
            
        };
    }
}

#endif /* ConstraintContact_h */
