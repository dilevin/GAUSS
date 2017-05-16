//
//  Scene.h
//  Gauss
//
//  Created by David Levin on 2/21/17.
//
//

#ifndef Scene_h
#define Scene_h

#include <World.h>

#include "PhysicsEntity.h"
#include "Update.h"

//FEM Renderables
#include "RenderableFEM.h"
#include "RenderableFEMHex.h"

#include <TimeStepper.h>

namespace Gauss {

    
    template<typename World, typename TimeStepper>
    class Scene {
        
    public:
        Scene(World *world) { std::cout<<"Wrong Scene Template \n"; }
    
    protected:
        
    private:
    };
    
    //QT Scene (template less base class which contains signals and slots
    class QtScene : public QObject {
        Q_OBJECT
    public:
        QtScene() { }
        virtual void initScene() = 0;
        virtual Qt3DCore::QEntity * getEntity() = 0;
        virtual void step() = 0;
        
    protected:
        
    private:
        
    public slots:
        void takeStep() { step(); }
    };
    
    template<typename DataType, typename ...SystemTypes,  typename TimeStepper>
    class Scene<World<DataType, SystemTypes...>, TimeStepper> : public QtScene {
        
    public:
        Scene(World<DataType,SystemTypes...> *world, TimeStepper *stepper) : QtScene() {
            m_world = world;
            m_stepper = stepper;
            m_rootEntity = NULL;
            initScene();
        }

        void initScene();
        
        Qt3DCore::QEntity * getEntity() {
            return m_rootEntity;
        }
        
        void step() {
            m_stepper->step(*m_world);
            
            for(auto &obj : m_entities) {
                obj->update(m_world->getState());
            }
        }
        
    protected:
        
        World<DataType, SystemTypes...> *m_world;
        std::vector<PhysicalEntity<DataType> *> m_entities;
        TimeStepper *m_stepper;
        Qt3DCore::QEntity *m_rootEntity;
        
    private:
        
    };
    
    template<typename DataType, typename ...SystemTypes, typename TimeStepper>
    void Scene<World<DataType, SystemTypes...>, TimeStepper >::initScene() {
        std::cout<<"Init Scene \n";
        
        if(!m_rootEntity)
            m_rootEntity = new Qt3DCore::QEntity();
              
        std::vector<PhysicalEntity<DataType> *> &list = m_entities;
        forEach(m_world->getSystemList(), [&list](auto a) {
            list.push_back(new PhysicalEntity<typename std::remove_pointer<decltype(a) >::type >(a) );
        });
    
        forEach(m_world->getForceList(), [&list](auto a) {
            list.push_back(new PhysicalEntity<typename std::remove_pointer<decltype(a) >::type >(a) );
        });
    
        
        //setup entities
        for(auto &entity : m_entities) {
            entity->getEntity(m_rootEntity, m_world->getState());
        }
        
    }
    
}


#endif /* Scene_h */
