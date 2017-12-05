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

//Big thing to figure out next is how to distinguish between getting DOFS supporting an vertex or those supporting an element
namespace Gauss {
    namespace Collisions {
        template<typename DataType, typename World, typename CollisionDetectorImpl>
        class CollisionDetector
        {
        public:
            
            template <typename ...Params>
            CollisionDetector(World &world, Params ...params) : m_impl(world, params...), m_world(world) {
            
            }
            
            ~CollisionDetector() { }
            
            //run this to get ready
            inline void detectCollisions() {
                m_impl.detectCollisions(m_world);
            }
            
            inline unsigned int getNumCollisions() {
                return m_impl.getNumCollisions();
            }
            
            //get Gradients
            template<typename Assembler, unsigned int Operation>
            inline void getGradient(Assembler &assembler, const State<DataType> &state) {
                
                ConstraintIndex index(0,0,1);
                
                auto &objAList = m_impl.getCollisionsObjectA();
                auto &objBList = m_impl.getCollisionsObjectB();
                auto &sharedInfo = m_impl.getSharedInfo();
                
                World &world =  m_world;
                
                //forEach loop for multivector
                forEach(objAList, [&assembler, &sharedInfo, &index,&world](auto &collisionInfo) {
                    std::cout<<"\n"<<sharedInfo[collisionInfo.getShared()].getNormal()<<"\n";
                    apply(world.getSystemList(), collisionInfo.getObject(), [&assembler, &collisionInfo, &index, &world, &sharedInfo](auto &a) {
                        auto J = sharedInfo[collisionInfo.getShared()].getNormal().transpose()*a->getDVDQ(sharedInfo[collisionInfo.getShared()].getPosition(), collisionInfo.getData(collisionInfo.collisionType));
                        assign(assembler, J, std::array<ConstraintIndex,1>{{index}},a->getQDot(collisionInfo.getData(collisionInfo.collisionType)));
                        });
                        
                        index.offsetGlobalId(1);
                
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
            World &m_world;
            
            
        private:
        };
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
