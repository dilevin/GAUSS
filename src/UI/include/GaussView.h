//
//  GaussView.h
//  Gauss
//
//  Created by David Levin on 4/8/17.
//
//

#ifndef GaussView_h
#define GaussView_h

#include <Qt3DIncludes.h>
#include <Scene.h>

#define GAUSSVIEW(scene)  Qt3DExtras::Qt3DWindow view; \
                          GaussView<std::remove_pointer<decltype(scene)>::type> gaussView(view, scene); \
                          view.show(); 

namespace Gauss {
    
        //wrapper class to deal with setup code for rendering scenes
        template<typename SceneType>
        class GaussView
        {
        public:
            GaussView(Qt3DExtras::Qt3DWindow  &view, SceneType *scene) {
                
                m_scene=scene;
                m_scene->takeStep(); //fix for the fact that my initialization doesn't seem to be kicking in.
                Qt3DCore::QEntity *sceneEntity = m_scene->getEntity();
                
                // Camera
                m_camera = view.camera();
                m_camera->lens()->setPerspectiveProjection(10.0f, 16.0f/9.0f, 0.1, 1000.0f);
                m_camera->setPosition(QVector3D(0, 0, 5.0f));
                m_camera->setViewCenter(QVector3D(0, 0, 0));
                
                // For camera controls
                m_camController = new Qt3DExtras::QOrbitCameraController(sceneEntity);
                m_camController->setLinearSpeed( 50.0f );
                m_camController->setLookSpeed( 180.0f );
                m_camController->setCamera(m_camera);
                
                view.setRootEntity(sceneEntity);
                //view.show();
                
                //Simple Timer to Play the simulation
                m_timer = new QTimer();
                QObject::connect(m_timer, SIGNAL(timeout()), scene, SLOT(takeStep()));
                
            }
            
            ~GaussView() {
                delete m_scene;
            }
            
            inline void startScene() { m_timer->start(); }
            inline void stopScene() { m_timer->stop(); }
            
        protected:
            
            SceneType *m_scene;
            Qt3DRender::QCamera *m_camera;
            Qt3DExtras::QOrbitCameraController *m_camController;
            QTimer *m_timer;
            
        private:
        };
}
#endif /* GaussView_h */
