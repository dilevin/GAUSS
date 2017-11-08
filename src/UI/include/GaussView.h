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
#include <QKeyEvent>
#include <Scene.h>

#define GAUSSVIEW(scene)  Gauss3DWindow view; \
                          GaussView<std::remove_pointer<decltype(scene)>::type> gaussView(view, scene); \
                          view.show(); 
class Gauss3DWindow: public Qt3DExtras::Qt3DWindow
{
    Q_OBJECT
public:
    Gauss3DWindow(QScreen *screen=nullptr):Qt3DExtras::Qt3DWindow(screen)
    {
        
    }
    
    ~Gauss3DWindow()
    {
        
    }
    
    void setSceneTimer(QTimer *sceneTimer) {
        m_timer = sceneTimer;
    }
    
    void startScene() { m_timer->start(); }
    void stopScene() { m_timer->stop(); }
    
protected:
    
    
    QTimer *m_timer;
    
    void keyPressEvent(QKeyEvent *ev)
    {
        switch (ev->key()) {
            case Qt::Key_P:
                m_timer->start();
                break;
            default:
                break;
        }
    }
};

namespace Gauss {
    
    
    
        //wrapper class to deal with setup code for rendering scenes
        template<typename SceneType>
        class GaussView
        {
        public:
            GaussView(Gauss3DWindow  &view, SceneType *scene) : m_view(view) {
                
                m_scene=scene;
                m_scene->takeStep(); //fix for the fact that my initialization doesn't seem to be kicking in.
                Qt3DCore::QEntity *sceneEntity = m_scene->getEntity();
                
                // Camera
                m_camera = m_view.camera();
                m_camera->lens()->setPerspectiveProjection(10.0f, 16.0f/9.0f, 0.1, 1000.0f);
                m_camera->setPosition(QVector3D(0, 0, 5.0f));
                m_camera->setViewCenter(QVector3D(0, 0, 0));
                
                // For camera controls
                m_camController = new Qt3DExtras::QOrbitCameraController(sceneEntity);
                m_camController->setLinearSpeed( 50.0f );
                m_camController->setLookSpeed( 180.0f );
                m_camController->setCamera(m_camera);
                
                m_view.setRootEntity(sceneEntity);
                
                //view.show();
                
                //Simple Timer to Play the simulation
                QTimer *timer = new QTimer();
                QObject::connect(timer, SIGNAL(timeout()), scene, SLOT(takeStep()));
                m_view.setSceneTimer(timer);
                
                
            }
            
            ~GaussView() {
                delete m_scene;
            }
            
            
            inline void startScene() {
                m_view.startScene();
            }
            
            inline void stopScene() {
                m_view.stopScene();
            }
            
        protected:
            
            SceneType *m_scene;
            Qt3DRender::QCamera *m_camera;
            Qt3DExtras::QOrbitCameraController *m_camController;
            Gauss3DWindow &m_view;

            
        private:
        };
}
#endif /* GaussView_h */
