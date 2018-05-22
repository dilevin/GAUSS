//
//  RenderableRigidBody.h
//  Gauss
//
//  Created by David Levin on 5/21/18.
//

#ifndef RenderableRigidBody_h
#define RenderableRigidBody_h

//
//  RenderableEmbeddedSurface.h
//  Gauss
//
//  Created by David Levin on 3/16/18.
//
#include <PhysicalSystemRigidBody.h>

namespace Gauss {
    //Renderable for a hexahedral finite element mesh
    template<typename DataType>
    class Renderable< RigidBodies::PhysicalSystemRigidBody<DataType> >: public Renderable<DataType>
    {
    public:
        
        
        Renderable(RigidBodies::PhysicalSystemRigidBody<DataType> *mesh) : Renderable<DataType>() {
            std::cout<<"Renderable PhysicalSystemRigidBody \n";
            m_mesh = mesh;
        }
        
        Qt3DCore::QEntity * getEntity(Qt3DCore::QEntity *root, State<DataType> &state) {
            
            //custom wire frame material
            Qt3DRender::QMaterial *wireframeMaterial = new Qt3DRender::QMaterial;
            Qt3DRender::QEffect *wireframeEffect = new Qt3DRender::QEffect;
            
            //parameters
            Qt3DRender::QParameter *m_ambient = new Qt3DRender::QParameter(QStringLiteral("ka"), QColor::fromRgbF(0.45f, 0.45f, 0.45f, 1.0f));
            Qt3DRender::QParameter *m_diffuse = new Qt3DRender::QParameter(QStringLiteral("kd"), QColor::fromRgbF(1.0f, 1.0f, 1.0f, 0.0f));
            Qt3DRender::QParameter *m_specular = new Qt3DRender::QParameter(QStringLiteral("ks"), QColor::fromRgbF(0.1f, 0.1f, 0.1f, 1.0f));
            Qt3DRender::QParameter *m_shininess = new Qt3DRender::QParameter(QStringLiteral("shininess"), 10.0);
            Qt3DRender::QParameter *m_lineWidth = new Qt3DRender::QParameter(QStringLiteral("line.width"), 3.0);
            Qt3DRender::QParameter *m_lineColor = new Qt3DRender::QParameter(QStringLiteral("line.color"), QColor::fromRgbF(0.0f, 0.0f, 0.0f, 1.0f));
            
            Qt3DRender::QTechnique *wfGL3Technique = new Qt3DRender::QTechnique;
            Qt3DRender::QRenderPass *wfGL3RenderPass = new Qt3DRender::QRenderPass;
            Qt3DRender::QShaderProgram *wfGL3Shader = new Qt3DRender::QShaderProgram;
            Qt3DRender::QFilterKey *filterKey = new Qt3DRender::QFilterKey;
            
            //load shader code
            wfGL3Shader->setVertexShaderCode(Qt3DRender::QShaderProgram::loadSource(QUrl::fromLocalFile(QString::fromStdString(dataDir())+"/shaders/robustwireframe.vert")));
            
            wfGL3Shader->setFragmentShaderCode(Qt3DRender::QShaderProgram::loadSource(QUrl::fromLocalFile(QString::fromStdString(dataDir())+"/shaders/robustwireframe.frag")));
            
            wfGL3Shader->setGeometryShaderCode(Qt3DRender::QShaderProgram::loadSource(QUrl::fromLocalFile(QString::fromStdString(dataDir())+"/shaders/robustwireframe.geom")));
            
            //filter decides when to use each shader
            filterKey->setName(QStringLiteral("renderingStyle"));
            filterKey->setValue(QStringLiteral("forward"));
            
            //setup render pass
            wfGL3RenderPass->setShaderProgram(wfGL3Shader);
            
            //setup technique
            wfGL3Technique->addFilterKey(filterKey);
            wfGL3Technique->addRenderPass(wfGL3RenderPass);
            wfGL3Technique->graphicsApiFilter()->setMajorVersion(3);
            wfGL3Technique->graphicsApiFilter()->setMinorVersion(0);
            wfGL3Technique->graphicsApiFilter()->setProfile(Qt3DRender::QGraphicsApiFilter::CoreProfile);
            
            //setup effect
            wireframeEffect->addParameter(m_ambient);
            wireframeEffect->addParameter(m_diffuse);
            wireframeEffect->addParameter(m_specular);
            wireframeEffect->addParameter(m_shininess);
            wireframeEffect->addParameter(m_lineWidth);
            wireframeEffect->addParameter(m_lineColor);
            wireframeEffect->addTechnique(wfGL3Technique);
            
            wireframeMaterial->setEffect(wireframeEffect);
            
            //############### END MATERIAL SETUP #################
            //custom geometry which stores geo/texture/normal/color data, data layout and property names
            //custom mesh renderer takes in geometry as well as information about offset per primitive and type of primitive to render
            m_meshRenderer = new Qt3DRender::QGeometryRenderer;
            
            //custom mesh renderer takes in geometry as well as information about offset per primitive and type of primitive to render
            m_meshGeometry = new Qt3DRender::QGeometry(m_meshRenderer);
            
            m_vertexDataBuffer = new Qt3DRender::QBuffer(Qt3DRender::QBuffer::VertexBuffer, m_meshGeometry);
            
            Qt3DRender::QBuffer *indexDataBuffer = new Qt3DRender::QBuffer(Qt3DRender::QBuffer::IndexBuffer, m_meshGeometry);
            
            //total numebr of verts = numVetsPerElement*numElements
            
            long totalScalarV = 3*m_mesh->getGeometry().second.rows()*m_mesh->getGeometry().second.cols();
            
            
            long numElements = m_mesh->getGeometry().second.rows();
            
            std::cout<<"Num Verts: "<<totalScalarV<<" Num Faces: "<<numElements<<"\n";
            
            //data for vertex positions + vertex normals + vertex colors = 3*totalScalarV
            m_vertexBufferData.resize(3*totalScalarV* sizeof(float)); //array for drawing
            m_updateData.resize(totalScalarV*sizeof(float)); //array for updating during animation
            
            
            //per face color and per face normals
            
            // Colors
            QVector3D red(1.0f, 1.0f, 0.0f);
            
            float *rawVertexArray = reinterpret_cast<float *>(m_vertexBufferData.data());
            float *rawUpdatePositionArray = reinterpret_cast<float *>(m_updateData.data());
            
            int idx = 0;
            int vertexId = 0;
            int elId = 0;
            const auto & V = m_mesh->getGeometry().first;
            const auto & F = m_mesh->getGeometry().second;
            
            
            //vertex positions
            for(; idx<totalScalarV; ++elId) {
                
                for(unsigned int ii=0; ii<F.cols(); ++ii) {
                    vertexId = F(elId, ii);
                    rawUpdatePositionArray[idx] = V(vertexId,0);
                    rawVertexArray[idx++] = V(vertexId,0);
                    rawUpdatePositionArray[idx] = V(vertexId,1);
                    rawVertexArray[idx++] = V(vertexId,1);
                    rawUpdatePositionArray[idx] = V(vertexId,2);
                    rawVertexArray[idx++] = V(vertexId,2);
                }
                
            }
            
            QVector3D v0, v1, v2;
            QVector3D n012;
            
            // Vector Normals
            QVector3D n[3];
            
            
            //really lazy vertex normals for now
            elId = 0;
            
            for(; elId < numElements; ) {
                
                v0 = QVector3D(V(F(elId,0),0), V(F(elId,0),1), V(F(elId,0),2));
                v1 = QVector3D(V(F(elId,1),0), V(F(elId,1),1), V(F(elId,1),2));
                v2 = QVector3D(V(F(elId,2),0), V(F(elId,2),1), V(F(elId,2),2));
                
                //0,1,2
                n012 = QVector3D::normal(v0, v1, v2);
                
                // Vector Normals
                n[0] = n012;
                n[1] = n012;
                n[2] = n012;
                
                for(unsigned int ii=0; ii< F.cols(); ++ii) {
                    rawVertexArray[idx++] = n[ii].x();
                    rawVertexArray[idx++] = n[ii].y();
                    rawVertexArray[idx++] = n[ii].z();
                }
                
                elId++;
            }
            
            elId = 0;
            
            //face colors
            for(; elId < numElements; ++elId) {
                
                for(unsigned int ii=0; ii<F.cols(); ++ii) {
                    rawVertexArray[idx++] = red.x();
                    rawVertexArray[idx++] = red.y();
                    rawVertexArray[idx++] = red.z();
                }
                
            }
            
            
            
            QByteArray indexBufferData;
            indexBufferData.resize(3*numElements*sizeof(ushort)); //ah indices, for sweet sweet index rendering
            ushort *rawIndexArray = reinterpret_cast<ushort *>(indexBufferData.data());
            
            idx = 0;
            elId = 0;
            
            for(unsigned int idy = 0;idx < 3*numElements; ) {
                
                rawIndexArray[idx++] = idy+0;
                rawIndexArray[idx++] = idy+1;
                rawIndexArray[idx++] = idy+2;
                
                idy += 3;
            }
            
            m_vertexDataBuffer->setData(m_vertexBufferData);
            indexDataBuffer->setData(indexBufferData);
            
            //entity takes in mesh renderer and uses it to draw the mesh
            // Attributes
            Qt3DRender::QAttribute *positionAttribute = new Qt3DRender::QAttribute();
            positionAttribute->setAttributeType(Qt3DRender::QAttribute::VertexAttribute);
            positionAttribute->setBuffer(m_vertexDataBuffer);
            positionAttribute->setVertexBaseType(Qt3DRender::QAttribute::Float);
            positionAttribute->setVertexSize(3);
            positionAttribute->setByteOffset(0);
            positionAttribute->setByteStride(0);
            positionAttribute->setCount(totalScalarV/3);
            positionAttribute->setName(Qt3DRender::QAttribute::defaultPositionAttributeName());
            
            Qt3DRender::QAttribute *normalAttribute = new Qt3DRender::QAttribute();
            normalAttribute->setAttributeType(Qt3DRender::QAttribute::VertexAttribute);
            normalAttribute->setBuffer(m_vertexDataBuffer);
            normalAttribute->setVertexBaseType(Qt3DRender::QAttribute::Float);
            normalAttribute->setVertexSize(3);
            normalAttribute->setByteOffset(totalScalarV * sizeof(float));
            normalAttribute->setByteStride(0);
            normalAttribute->setCount(totalScalarV/3);
            normalAttribute->setName(Qt3DRender::QAttribute::defaultNormalAttributeName());
            
            Qt3DRender::QAttribute *colorAttribute = new Qt3DRender::QAttribute();
            colorAttribute->setAttributeType(Qt3DRender::QAttribute::VertexAttribute);
            colorAttribute->setBuffer(m_vertexDataBuffer);
            colorAttribute->setVertexBaseType(Qt3DRender::QAttribute::Float);
            colorAttribute->setVertexSize(3);
            colorAttribute->setByteOffset((2*totalScalarV)* sizeof(float));
            colorAttribute->setByteStride(0);
            colorAttribute->setCount(totalScalarV/3);
            colorAttribute->setName(Qt3DRender::QAttribute::defaultColorAttributeName());
            
            Qt3DRender::QAttribute *indexAttribute = new Qt3DRender::QAttribute();
            indexAttribute->setAttributeType(Qt3DRender::QAttribute::IndexAttribute);
            indexAttribute->setBuffer(indexDataBuffer);
            indexAttribute->setVertexBaseType(Qt3DRender::QAttribute::UnsignedShort);
            indexAttribute->setVertexSize(1);
            indexAttribute->setByteOffset(0);
            indexAttribute->setByteStride(0);
            indexAttribute->setCount(totalScalarV/3);
            
            m_meshGeometry->addAttribute(positionAttribute);
            m_meshGeometry->addAttribute(normalAttribute);
            m_meshGeometry->addAttribute(colorAttribute);
            m_meshGeometry->addAttribute(indexAttribute);
            
            m_meshRenderer->setInstanceCount(1);
            m_meshRenderer->setIndexOffset(0);
            m_meshRenderer->setFirstInstance(0);
            m_meshRenderer->setPrimitiveType(Qt3DRender::QGeometryRenderer::Triangles);
            m_meshRenderer->setGeometry(m_meshGeometry);
            
            m_meshRenderer->setVertexCount(3*numElements);
            
            //transform
            m_transform = new Qt3DCore::QTransform();
            m_transform->setScale3D(QVector3D(1,1,1));
            m_transform->setRotationX(0.0);
            m_transform->setRotationY(0.0);
            m_transform->setRotationZ(0.0);
            m_transform->setTranslation(QVector3D(0,0,0));
            
            //Qt3DExtras::QPhongMaterial *material = new Qt3DExtras::QPhongMaterial();
            //Qt3DRender::QMaterial *material = new Qt3DExtras::QPerVertexColorMaterial();
            
            //material->setDiffuse(QColor(220.0, 0.0,0.0));*/
            
            m_entity = new Qt3DCore::QEntity(root);
            m_entity->addComponent(m_meshRenderer);
            m_entity->addComponent(m_transform);
            m_entity->addComponent(wireframeMaterial);
            
            
            
            return m_entity;
            
        }
        
        void update(State<DataType> &state) const {
            
            //reset mesh positions and add in the state (which I think is already setup in the right order ?)
            
            
            //1.) Grab raw data
            float *rawData = const_cast<float *>(reinterpret_cast<const float *>(m_updateData.data()));
            //float *rawOriginalData = const_cast<float *>(reinterpret_cast<const float *>(m_vertexBufferData.data()));
            
            //2.) Update using mesh data + state which for linear FEM is the displacement of the node
            //Eigen::Map<Eigen::VectorXd> pos = mapDOFEigen(m_fem->getQ(), state);
            Eigen::Vector3x<DataType> pos;
            
            unsigned int vertexId, idx;
            idx = 0;
            
           for(unsigned int elId=0;  elId < m_mesh->getGeometry().second.rows(); ++elId) {
                
                for(unsigned int ii=0; ii<m_mesh->getGeometry().second.cols(); ++ii) {
                    
                    
                    vertexId = m_mesh->getGeometry().second(elId, ii);
                    pos = m_mesh->getPosition(state, vertexId);
                    rawData[idx] = pos[0];
                    idx++;
                    rawData[idx] = pos[1];
                    idx++;
                    rawData[idx] = pos[2];
                    idx++;
                }
                
            }
            
            //Update
            m_vertexDataBuffer->updateData(0, m_updateData);
            
        }
        
    protected:
        Qt3DCore::QEntity *m_entity;
        Qt3DRender::QGeometryRenderer *m_meshRenderer;
        Qt3DRender::QGeometry *m_meshGeometry;
        RigidBodies::PhysicalSystemRigidBody<DataType> *m_mesh;
        Qt3DCore::QTransform *m_transform;
        QByteArray m_updateData;
        QByteArray m_vertexBufferData;
        Qt3DRender::QBuffer *m_vertexDataBuffer;
        
    private:
    };
    
}

#endif /* RenderableRigidBody_h */
