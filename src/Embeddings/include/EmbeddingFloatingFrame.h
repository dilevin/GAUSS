//
//  EmbeddingFloatingFrame.h
//  Gauss
//
//  Created by David Levin on 5/24/18.
//

#ifndef EmbeddingFloatingFrame_h
#define EmbeddingFloatingFrame_h

#include <GaussIncludes.h>
#include <FEMIncludes.h>
#include <PhysicalSystemRigidBody.h>

namespace Gauss {
    namespace Embeddings {
        template<typename DataType, typename System>
        class  PhysicalSystemFloatingFrameImpl : public PhysicalSystemRigidBodyImpl<DataType> {
            
        public:
            
            
            using PhysicalSystemImpl = typename System::ImplType;
            using Embedding = EmbeddingFunction<DataType,PhysicalSystemImpl>;
            using Embedding::PhysicalSystemImpl::getEnergy;
            using Embedding::PhysicalSystemImpl::getStrainEnergy;
            using Embedding::PhysicalSystemImpl::getStrainEnergyPerElement;
            using Embedding::PhysicalSystemImpl::getMassMatrix;
            using Embedding::PhysicalSystemImpl::getStiffnessMatrix;
            using Embedding::PhysicalSystemImpl::getForce;
            using Embedding::PhysicalSystemImpl::getInternalForce;
            using Embedding::PhysicalSystemImpl::getQ;
            using Embedding::PhysicalSystemImpl::getQDot;
            
            //contructor
            PhysicalSystemFloatingFrameImpl(PhysicalSystem<DataType, System> *embeded) : PhysicalSystemRigidBodyImpl(params...) {
                
                //initialize Rigid body from V,F of Physical System, I'm going to assume it's not an embedded mesh and initalie off of the boundary facets
                
                
                //add this and fem to a world
            }
            
            //Geometry stuff
            /*template<typename ...Params>
            inline decltype(auto) getPosition(const State<DataType> &state, Params &...params) const {
                return m_embedding.getPosition(*this, state, params...);
            }
            
            template<typename ...Params>
            inline decltype(auto) getDPDQ(const State<DataType> &state, Params &...params) {
                return m_embedding.getDPDQ(*this, state, params...);
            }
            
            template<typename ...Params>
            inline decltype(auto) getVelocity(const State<DataType> &state, Params &...params) {
                return m_embedding.getVelocity(*this, state, params...);
            }
            
            template<typename ...Params>
            inline decltype(auto) getDVDQ(const State<DataType> &state, Params &...params) {
                return m_embedding.getDVDQ(*this, state, params...);
            }
            
            //Get function supporing point in space
            template<typename Vector>
            inline  decltype(auto) getQ(const Vector &x, unsigned int elementId) const {
                return m_embedding.getQ(*this, x, elementId);
            }
            
            //Get function supporing vertex in space
            inline  decltype(auto) getQ(unsigned int elementId) const {
                return m_embedding.getQ(*this, elementId);
            }
            
            //Get function supporing vertex in space
            inline  decltype(auto) getQDot(unsigned int elementId)  const {
                
                return m_embedding.getQDot(*this, elementId);
            }
            
            inline decltype(auto) getGeometry() { return m_embedding.getGeometry(); }*/
            
        protected:
            
            DOFPair<decltype(getQ()) &, decltype(m_embeded->getQ()) &> m_q;
            DOFPair<decltype(getQDot()) &, decltype(m_embeded->getQDot()) &> m_qDot;
            PhysicalSystem<DataType, System> *m_embeded
            
        private:
            
            
        };
        
        template<typename DataType, typename System>
        using PhysicalSystemFloatingFrame = PhysicalSystem<DataType, PhysicalSystemFloatingFrameImpl<DataType, System> >;
    }
}
#endif /* EmbeddingFloatingFrame_h */
