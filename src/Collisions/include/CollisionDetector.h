//
//  CollisionDetector.h
//  Collision detector parent class
//
//  Created by David Levin on 4/20/17.
//
//

#ifndef CollisionDetector_h
#define CollisionDetector_h

#include <Constraint.h>
#include <UtilitiesEigen.h>
#include <GaussIncludes.h>

//This should become a constraint (just makes everything easier)
namespace Gauss {
    namespace Collisions {
        template<typename DataType, typename CollisionDetectorImpl>
        class CollisionDetector
        {
        public:
            
            template <typename World, typename ...Params>
            CollisionDetector(World &world, Params ...params) : m_impl(world, params...) {
            
            }
            
            ~CollisionDetector() { }
            
            //run this to get ready
            //if this is a constraint how do I call detect collisions
            //send extra reference to any thing that needs access to the collision detector
            //collision detectors are in constraint list and collision detection list
            //
            template<typename World>
            inline void update(World &world) {
                detectCollisions(world);
            }
            
            template<typename World>
            inline void detectCollisions(World &world) {
                m_impl.detectCollisions(world);
            }
            
            inline unsigned int getNumRows() {
                return getNumCollisions();
            }
            
            inline unsigned int getNumCollisions() {
                return m_impl.getNumCollisions();
            }
            
            //get Gradients
            template<typename World, typename Assembler, unsigned int Operation>
            void getGradient(Assembler &assembler, World &world, const State<DataType> &state, const ConstraintIndex &index) {
                
                ConstraintIndex indexInc = index;
                indexInc.setNumRows(1);
                
                auto &objAList = m_impl.getCollisionsObjectA();
                auto &objBList = m_impl.getCollisionsObjectB();
                auto &sharedInfo = m_impl.getSharedInfo();
                
                //forEach loop for multivector
                forEach(objAList, [&assembler, &sharedInfo, &indexInc,&world, &state](auto &collisionInfo) {
                    apply(world.getSystemList(), collisionInfo.getObject(), [&assembler, &collisionInfo, &indexInc, &world, &sharedInfo, &state](auto &a) {
                        
                        //Build Jacobian differently for vertex collisions and in element collision events
                        static_if<std::remove_reference<decltype(collisionInfo)>::type::collisionType == 0>([&](auto f){
                            auto J = (sharedInfo[collisionInfo.getShared()].getNormal().transpose()*a->getDPDQ(state, collisionInfo.getData(collisionInfo.collisionType))).eval();
                            assign(assembler, J, std::array<ConstraintIndex,1>{{indexInc}},a->getQDot(collisionInfo.getData(collisionInfo.collisionType)));
                        }).else_([&](auto f) {
                        
                            auto J = (sharedInfo[collisionInfo.getShared()].getNormal().transpose()*a->getDPDQ(state, collisionInfo.getData(collisionInfo.collisionType),sharedInfo[collisionInfo.getShared()].getPosition())).eval();
                            auto q = a->getQDot(sharedInfo[collisionInfo.getShared()].getPosition(), collisionInfo.getData(collisionInfo.collisionType));
                            //std::cout<<"Q SIZE: "<<q.size()<<"\n";
                            assign(assembler, J, std::array<ConstraintIndex,1>{{indexInc}},q);
                        });
                        
                        
                    });
                    
                    //need another loop of objectB
                    indexInc.offsetGlobalId(1);
                
                });
                
                indexInc = index;
                indexInc.setNumRows(1);
                
                //forEach loop for multivector
                forEach(objBList, [&assembler, &sharedInfo, &indexInc,&world, &state](auto &collisionInfo) {
                    if(collisionInfo.getObject().type() != -1) {
                        apply(world.getSystemList(), collisionInfo.getObject(), [&assembler, &collisionInfo, &indexInc, &world, &sharedInfo, &state](auto &a) {
                            
                            //Build Jacobian differently for vertex collisions and in element collision events
                            static_if<std::remove_reference<decltype(collisionInfo)>::type::collisionType == 0>([&](auto f){
                                auto J = (-sharedInfo[collisionInfo.getShared()].getNormal().transpose()*a->getDPDQ(state, collisionInfo.getData(collisionInfo.collisionType))).eval();
                                assign(assembler, J, std::array<ConstraintIndex,1>{{indexInc}},a->getQ(collisionInfo.getData(collisionInfo.collisionType)));
                            }).else_([&](auto f) {
                                
                                //something is wrong here
                                auto J = (-sharedInfo[collisionInfo.getShared()].getNormal().transpose()*a->getDPDQ(state, collisionInfo.getData(collisionInfo.collisionType),sharedInfo[collisionInfo.getShared()].getPosition())).eval();
                                auto q = a->getQ(sharedInfo[collisionInfo.getShared()].getPosition(), collisionInfo.getData(collisionInfo.collisionType));
                                //std::cout<<"Q SIZE: "<<q.size()<<"\n";
                                assign(assembler, J, std::array<ConstraintIndex,1>{{indexInc}},q);
                            });
                            
                            
                        });
                    }
                    //need another loop of objectB
                    indexInc.offsetGlobalId(1);
                    
                });
                
                //forLoop<IsParallel<Assembler>::value>(objAList, assembler, [&world, &index, &cd](auto &assemble, auto &constraint) {
                    
                    //TODO add static if inside apply to handle different constraint types.
                    //apply(world, constraint.getObject(), [&assemble, &constraint, &index](auto a) {
                      //  auto J = constraint.getShared().getNormal().Transpose()*a->getDVDQ(constraint.getData(1));
                        //assign(assemble, J, std::array<ConstraintIndex,1>{{index}},a->getQDotDOF(constraint.getData(1)));
                    //});
                    
                    //index.offsetGlobalId(1);
                    //std::cout<<constraint.getPositon()<<"\n";
                    
                //});
                
                //handle object B collisions
                //for each loop for multivector
            }
        
            //get Values
            template<typename Assembler>
            inline void getFunction(Assembler &assembler, const State<DataType> &state) const {
                std::cout<<"Not implemented yet \n";
                assert(0 ==1);
            }
            
        protected:
            
            CollisionDetectorImpl m_impl;
        
            
        private:
        };
        
        template<typename DataType, template<typename A> class DetectorImpl>
        using ConstraintCollisionDetector = Constraint<DataType, CollisionDetector<DataType, DetectorImpl<DataType> > >;
    }
}

/*
 * contact constraint
 * (pa(q(t+dt)) - pb(q(t+dt)))*n > 0
 * constraint Jacobian [n'*dpa/dq,  -n'*dpb/dq] [dqa/dt dtdqb/dt*dt]'
 * so I need methods that can compute dpa/dq given position on an objet
 * leaning towards constraint detector becoming inequality constraint again, get gradient methods should take in a constraintIndex
 */
#endif /* ConstraintCollisions_h */
