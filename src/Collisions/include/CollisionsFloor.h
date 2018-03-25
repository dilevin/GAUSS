//Just making some notes for myself
//template<BLORP>
/*planeDistance(object, plane) {
 doNothing
 */
//Requirements of object to work with floor collision
//Needs a method that tells you how far way it is from a plane
//planeDist(obj, plane params)
//Different planeDists for different geometry representations
//Need method to check geometry representation of an object
//Can use trait classes to check for certain methods in the implementation
// for i = 1 to numObjects
//      planeDistance(Object)
// end
// collision constraint needs to compute dP/dq|(collision position)
// how I access collision position depends on info I have
// have different types of collision constraints
// ideally those dont't live here, they live in the parent collision detector class

//compile time if
//if<getV>
//planeDist(system, plane)
//else geometry
//planeDist
//else fail
#ifndef _COLLISIONFLOOR_H
#define _COLLISIONFLOOR_H

#include <GaussIncludes.h>
#include <UtilitiesContact.h>

namespace Gauss {
    namespace Collisions {
        
        //I've made a lot of different attempts at getting this right (sgh)
        template<typename Geometry>
        class DetectCollisions {
        public:
            template<typename DataType, typename Object, typename ...Types>
            inline static void planeContact(int type, int index, const Object *obj, const State<DataType> &state,
                                       MultiVector<Types...> &listA, MultiVector<Types...> &listB,
                                       std::vector<SharedCollisionInfo<DataType> > &info,
                                       Eigen::Vector3x<DataType> &normal,
                                       Eigen::Vector3x<DataType> &pos) {
                std::cout<<"Generic planeContact routine does nothing \n";
                exit(0);
            }
            
        protected:
        private:
        };
        
        //Contact routines for vertex/element geometries
        template<>
        class DetectCollisions<std::pair<Eigen::MatrixXd &, Eigen::MatrixXi &> >{
        public:
            template<typename DataType, typename Object, typename ...Types>
            inline static void planeContact(int type, int index, const Object *obj, const State<DataType> &state,
                                       MultiVector<Types...> &listA, MultiVector<Types...> &listB,
                                       std::vector<SharedCollisionInfo<DataType> > &info,
                                       Eigen::Vector3x<DataType> &normal,
                                       Eigen::Vector3x<DataType> &pos) {
                
                for(unsigned int iv =0; iv < obj->getImpl().getV().rows(); ++iv) {
                    //if the object is on the wrong side of the floor, mark the collision
                    double proj = (obj->getPosition(state, iv) - pos).dot(normal);
                    
                    if(proj < 0){
                        //store collision
                        info.push_back(SharedCollisionInfo<DataType>(normal,obj->getPosition(state, iv)));
                        listA.add(ObjectCollisionInfo<DataType,0>(info.size()-1, SystemIndex(type,index), iv));
                        listB.add(ObjectCollisionInfo<DataType,0>(info.size()-1, SystemIndex(-1,0), 0));
                        
                    }
                }
                
                
            }
            
        protected:
        private:
        };
        
        //Some useful methods for plane-object collision checking
        template<typename DataType>
        class CollisionFloorImpl {

        public:
            
            template<typename World>
            CollisionFloorImpl(World &world, Eigen::Vector3x<DataType> floorNormal, Eigen::Vector3x<DataType> floorPosition) {
                m_floorNormal = floorNormal;
                m_floorPosition = floorPosition;
            }
            
            ~CollisionFloorImpl() { }
            
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
            
            Eigen::Vector3x<DataType> m_floorNormal;
            Eigen::Vector3x<DataType> m_floorPosition;
            
            //only vertex collisions for now
            MultiVector<ObjectCollisionInfo<DataType, 0> > m_objAList; //collision info for object A
            MultiVector<ObjectCollisionInfo<DataType, 0> > m_objBList; //collision info for object B
            std::vector<SharedCollisionInfo<DataType>> m_sharedList; //normals and world space positions;
        private:
            
        };

    }
}

template<typename DataType>
template<typename World>
void Gauss::Collisions::CollisionFloorImpl<DataType>::detectCollisions(World &world) {
    
    m_sharedList.clear();
    m_objAList.clear();
    m_objBList.clear();
    
    MultiVector<Gauss::Collisions::ObjectCollisionInfo<DataType, 0> > &objAList = m_objAList;
    MultiVector<Gauss::Collisions::ObjectCollisionInfo<DataType, 0> > &objBList = m_objBList;
    std::vector<Gauss::Collisions::SharedCollisionInfo<DataType> > &sharedList = m_sharedList;
    Eigen::Vector3x<DataType> &floorNormal = m_floorNormal;
    Eigen::Vector3x<DataType> &floorPosition = m_floorPosition;
    
    //Loop through every object, check if any points are on the wrong side of the floor, if so
    //record collision
    forEachIndex(world.getSystemList(), [&objAList, &objBList, &sharedList, &floorNormal, &floorPosition, &world](auto type, auto index, auto &a) {
        DetectCollisions<decltype(a->getGeometry())>::planeContact(type, index, a, world.getState(), objAList, objBList, sharedList, floorNormal, floorPosition);
    });
    

    
}

#endif
