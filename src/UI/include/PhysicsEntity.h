#ifndef PHYSICSENTITY_H
#define PHYSICSENTITY_H


#include <State.h>
#include <Qt3DCore/qentity.h>
#include <Update.h>
#include <RenderableFEM.h>
#include <RenderableFEMHex.h>
#include <RenderableFEMTri.h>
#include <RenderableEmbeddedSurface.h>

//specific physical system types
#include <PhysicalSystemParticles.h>
#include <PhysicalSystemFEM.h>
#include <ForceExternal.h>
#include <ForceSpring.h>

namespace Gauss {

    //Interface for physical entity
    template<typename DataType, typename ...Types>
    class PhysicalEntity
    {
    public:
        PhysicalEntity() {
            m_rootEntity = NULL;
        }
        
        Qt3DCore::QEntity * getEntity(Qt3DCore::QEntity *root, State<DataType> &state) {
            
            if(!m_rootEntity)
                m_rootEntity = new Qt3DCore::QEntity(root);
            else
                m_rootEntity->setParent(root);
            
            //concatanate all entities to a root node and return
            for(auto &obj : m_renderableObjects) {
                obj->getEntity(m_rootEntity, state);
            }
            
            return m_rootEntity;
        }
        
        void update(State<DataType> &state) {
            for(auto &obj : m_renderableObjects) {
                obj->update(state);
            }
        }
        
    protected:
        
        Qt3DCore::QEntity *m_rootEntity;
        std::vector<Renderable<DataType> *> m_renderableObjects;
        
    private:
        
    };
    
    //Do nothing entity to handle cases that haven't been implemented yet
    template<typename DataType, typename Impl>
    class PhysicalEntity<Force<DataType, Impl> >:
    public PhysicalEntity<DataType>
    
    {
    public:
        PhysicalEntity(Force<DataType, Impl> *system) :
        PhysicalEntity<DataType>()
        {
        }
        
        
        
    protected:
        
        
    private:
        
    };

    
    template<typename DataType, typename Position1, typename Position2>
    class PhysicalEntity<ParticleSystem::ForceSpring<DataType, Position1, Position2> >:
    public PhysicalEntity<DataType>
    
    {
    public:
        PhysicalEntity(ParticleSystem::ForceSpring<DataType, Position1, Position2> *system) :
        PhysicalEntity<DataType>()
        {
            //render connection between dofs
            //using Pos0 = typename std::remove_reference<decltype(system->getImpl())>::Position0Type;
            //using Pos1 = typename std::remove_reference<decltype(system->getImpl())>::Position1Type;
            
            PhysicalEntity<DataType>::m_renderableObjects.push_back(
               new RenderableLine<DataType, Position1, Position2>(system->getImpl().getPosition0(), system->getImpl().getPosition1())
            );
        }
        
        
        
    protected:
        
        
    private:
        
    };
    
    //Particle System
    //Interface for physical entity
    template<typename DataType>
    class PhysicalEntity<ParticleSystem::PhysicalSystemParticleSingle<DataType> > :
    public PhysicalEntity<DataType>
    
    {
    public:
        PhysicalEntity(ParticleSystem::PhysicalSystemParticleSingle<DataType> *system) :
        PhysicalEntity<DataType>()
        {
            PhysicalEntity<DataType>::m_renderableObjects.push_back(new Renderable<typename std::remove_reference<decltype(system->getQ())>::type >(system->getQ()));
        }
        
       

    protected:
        
        
    private:
        
    };
    
    //Finite Elements
    template<typename DataType,  unsigned int Num,
    template<typename Type, typename Energy> class QuadT,
    template<typename Type, typename Energy> class QuadU,
    template<typename Type, typename Func> class KE,
    template<typename Type, typename Func> class PE,
    template<typename Type, typename Func> class BF,
    template<typename Type> class SF,
    template<   typename A, unsigned int N,
    template<typename Type, typename Energy> class QuadratureRuleT,
    template<typename Type, typename Energy> class QuadratureRuleU,
    template<typename Type, typename Func> class KineticEnergy,
    template<typename Type, typename Func> class PotentialEnergy,
    template<typename Type, typename Func> class BodyForce,
    template<typename Type> class ShapeFunction >
    class ElementBase >
    class PhysicalEntity<PhysicalSystem<DataType, FEM::PhysicalSystemFEMImpl<DataType, ElementBase<DataType, Num, QuadT, QuadU, KE, PE, BF, SF> > > > :
    public PhysicalEntity<DataType>
    {
    public:
        using FEMSystem = PhysicalSystem<DataType, FEM::PhysicalSystemFEMImpl<DataType, ElementBase<DataType, Num, QuadT, QuadU, KE, PE, BF, SF> > >;

        
        PhysicalEntity(FEMSystem *system) : PhysicalEntity<DataType>()
        {
            std::cout<<"FEM System \n";
            PhysicalEntity<DataType>::m_renderableObjects.push_back(new Renderable<typename std::remove_pointer<decltype(system)>::type >(system));
        }
        
        
        
    protected:
        
        
    private:
        
    };

    //Embedded Mesh
    template<typename DataType, typename System>
    class PhysicalEntity<Embeddings::PhysicalSystemEmbeddedMesh<DataType, System> > : public PhysicalEntity<DataType>
    
    {
    public:
        PhysicalEntity(Embeddings::PhysicalSystemEmbeddedMesh<DataType, System> *system) :
        PhysicalEntity<DataType>()
        {
            PhysicalEntity<DataType>::m_renderableObjects.push_back(new Renderable<typename std::remove_pointer<decltype(system)>::type>(system));
        }
        
        
        
    protected:
        
        
    private:
        
    };
    
    
}

#endif
