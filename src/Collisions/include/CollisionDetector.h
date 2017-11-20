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

namespace Gauss {
    namespace Collisions {
        template<typename DataType, typename CollisionDetectorImpl>
        class CollisionDetector
        {
        public:
            
            template <typename World, typename ...Params>
            CollisionDetector(World &world, Params ...params) : m_world(world), m_impl(m_world, params) {
            
            }
            
            ~CollisionDetector() { }
            
            //run this to get ready
            inline void detectCollisions(const State<DataType> &state) {
                m_impl.detectCollisions(m_world);
             
                m_constraintList.resize(m_impl.getNumContacts());
                
                for(unsigned int ii=0;ii<m_impl.getNumContacts(); ++ii) {
                    m_constraintList[ii] = m_impl.getContact(ii);
                }
            }
            
            //get Gradients
            template<typename Assembler>
            inline void getGradient(Assembler &assembler, const State<DataType> &state) const {
                
                forLoop<IsParallel<Assembler>::value>(m_constraintList, assembler, [&](auto &assemble, auto &constraint) {
                    constraint->getGradient(assemble, state);
                });
            }
            
            //get Values
            template<typename Assembler>
            inline void getFunction(Assembler &assembler, const State<DataType> &state) const {
                forLoop<IsParallel<Assembler>::value>(m_constrainList, assembler, [&](auto &assemble, auto &constraint) {
                    constraint->getFunction(assemble, state);
                });
            }
            
            
        protected:
            
            CollisionDetectorImpl m_impl;
            World &m_world;
            
            //collision constraints list
            std::vector<double> m_constraintList;
            
            
        private:
        };
    }
}

#endif /* ConstraintCollisions_h */
