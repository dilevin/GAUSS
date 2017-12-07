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
                
                //update the constraint numberings in the world here
            }
            
            inline unsigned int getNumRows() {
                return getNumCollisions();
            }
            
            inline unsigned int getNumCollisions() {
                return m_impl.getNumCollisions();
            }
            
            //get Gradients
            template<typename World, typename Assembler, unsigned int Operation>
            inline void getGradient(Assembler &assembler, World &world, const State<DataType> &state, const ConstraintIndex &index) {
                
                ConstraintIndex indexInc = index;
                indexInc.setNumRows(1);
                
                auto &objAList = m_impl.getCollisionsObjectA();
                auto &objBList = m_impl.getCollisionsObjectB();
                auto &sharedInfo = m_impl.getSharedInfo();
                
                //forEach loop for multivector
                forEach(objAList, [&assembler, &sharedInfo, &indexInc,&world](auto &collisionInfo) {
                    std::cout<<"\n"<<sharedInfo[collisionInfo.getShared()].getNormal()<<"\n";
                    apply(world.getSystemList(), collisionInfo.getObject(), [&assembler, &collisionInfo, &indexInc, &world, &sharedInfo](auto &a) {
                        auto J = sharedInfo[collisionInfo.getShared()].getNormal().transpose()*a->getDVDQ(sharedInfo[collisionInfo.getShared()].getPosition(), collisionInfo.getData(collisionInfo.collisionType));
                        assign(assembler, J, std::array<ConstraintIndex,1>{{indexInc}},a->getQDot(collisionInfo.getData(collisionInfo.collisionType)));
                        });
                        
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
