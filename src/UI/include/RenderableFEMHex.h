//
//  RenderableFEMHex.h
//  Gauss
//
//  Created by David Levin on 5/11/17.
//
//

#ifndef RenderableFEMHex_h
#define RenderableFEMHex_h

#include <FEMIncludes.h>

namespace Gauss {
    
    //Renderable for a hexahedral finite element mesh 
    template<typename DataType, template <typename D, typename E> class QuadratureRule, template <typename F, typename G> class EnergyPotential>
    class Renderable<PhysicalSystem<DataType, FEM::PhysicalSystemFEMImpl<DataType,
    FEM::Element<DataType, 8, QuadratureRule,
    FEM::EnergyKineticNonLumped,
    EnergyPotential,
    FEM::BodyForceGravity,
    FEM::ShapeFunctionHexTrilinear> > > >: public Renderable<DataType>
    {
    public:
        
        using FEMSystem = PhysicalSystem<DataType, FEM::PhysicalSystemFEMImpl<DataType,
        FEM::Element<DataType, 8, QuadratureRule,
        FEM::EnergyKineticNonLumped,
        EnergyPotential,
        FEM::BodyForceGravity,
        FEM::ShapeFunctionHexTrilinear> > >;
        
        
        Renderable(FEMSystem *fem) : Renderable<DataType>() {
            std::cout<<"Renderable Trilinear Hexahedral FEM \n";
            m_fem = fem;
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
            
            
            //long totalScalarV = m_fem->getImpl().getV().rows()*m_fem->getImpl().getV().cols();
            
            //total numebr of verts = numVetsPerElement*numElements
            long totalScalarV = 3*m_fem->getImpl().getF().cols()*m_fem->getImpl().getF().rows();
            long numElements = m_fem->getImpl().getF().rows();
            
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
            const auto & V = m_fem->getImpl().getV();
            const auto & F = m_fem->getImpl().getF();
            
            
            //vertex positions
            /*for(; idx < totalScalarV; ++vertexId) {
             rawUpdatePositionArray[idx] = m_fem->getImpl().getV()(vertexId,0);
             rawVertexArray[idx++] = m_fem->getImpl().getV()(vertexId,0);
             
             rawUpdatePositionArray[idx] = m_fem->getImpl().getV()(vertexId,1);
             rawVertexArray[idx++] = m_fem->getImpl().getV()(vertexId,1);
             
             rawUpdatePositionArray[idx] = m_fem->getImpl().getV()(vertexId,2);
             rawVertexArray[idx++] = m_fem->getImpl().getV()(vertexId,2);
             
             }*/
            
            for(; idx<totalScalarV; ++elId) {
                
                for(unsigned int ii=0; ii<F.cols(); ++ii) {
                    vertexId = F(elId, ii);
                    rawUpdatePositionArray[idx] = m_fem->getImpl().getV()(vertexId,0);
                    rawVertexArray[idx++] = m_fem->getImpl().getV()(vertexId,0);
                    rawUpdatePositionArray[idx] = m_fem->getImpl().getV()(vertexId,1);
                    rawVertexArray[idx++] = m_fem->getImpl().getV()(vertexId,1);
                    rawUpdatePositionArray[idx] = m_fem->getImpl().getV()(vertexId,2);
                    rawVertexArray[idx++] = m_fem->getImpl().getV()(vertexId,2);
                }
                
                /*)//vertex 2
                vertexId = F(elId, 1);
                rawUpdatePositionArray[idx] = m_fem->getImpl().getV()(vertexId,0);
                rawVertexArray[idx++] = m_fem->getImpl().getV()(vertexId,0);
                rawUpdatePositionArray[idx] = m_fem->getImpl().getV()(vertexId,1);
                rawVertexArray[idx++] = m_fem->getImpl().getV()(vertexId,1);
                rawUpdatePositionArray[idx] = m_fem->getImpl().getV()(vertexId,2);
                rawVertexArray[idx++] = m_fem->getImpl().getV()(vertexId,2);
                
                //vertex 3
                vertexId = F(elId, 2);
                rawUpdatePositionArray[idx] = m_fem->getImpl().getV()(vertexId,0);
                rawVertexArray[idx++] = m_fem->getImpl().getV()(vertexId,0);
                rawUpdatePositionArray[idx] = m_fem->getImpl().getV()(vertexId,1);
                rawVertexArray[idx++] = m_fem->getImpl().getV()(vertexId,1);
                rawUpdatePositionArray[idx] = m_fem->getImpl().getV()(vertexId,2);
                rawVertexArray[idx++] = m_fem->getImpl().getV()(vertexId,2);
                
                //vertex 4
                vertexId = F(elId, 3);
                rawUpdatePositionArray[idx] = m_fem->getImpl().getV()(vertexId,0);
                rawVertexArray[idx++] = m_fem->getImpl().getV()(vertexId,0);
                rawUpdatePositionArray[idx] = m_fem->getImpl().getV()(vertexId,1);
                rawVertexArray[idx++] = m_fem->getImpl().getV()(vertexId,1);
                rawUpdatePositionArray[idx] = m_fem->getImpl().getV()(vertexId,2);
                rawVertexArray[idx++] = m_fem->getImpl().getV()(vertexId,2);*/
                
            }
            
            QVector3D v0, v1, v2, v3, v4, v5, v6, v7;
            QVector3D n023;
            QVector3D n012;
            QVector3D n765,n754,n037,n074,n621,n561,n326,n367,n104,n145;
            
            // Vector Normals
            QVector3D n[8];
            
            
            //really lazy vertex normals for now
            elId = 0;
            
            for(; elId < numElements; ) {
                
                //for hexes, split each element into 2 triangles for now so this adds 12 triangles
                //    4-----------5
                //    /|         /|
                //   7-|--------6 |
                //   | |        | |
                //   | 0 -------|-1    +y
                //   |/         |/     |
                //   3----------2      0-- +x
                //                    /
                //                   +z
            
                v0 = QVector3D(V(F(elId,0),0), V(F(elId,0),1), V(F(elId,0),2));
                v1 = QVector3D(V(F(elId,1),0), V(F(elId,1),1), V(F(elId,1),2));
                v2 = QVector3D(V(F(elId,2),0), V(F(elId,2),1), V(F(elId,2),2));
                v3 = QVector3D(V(F(elId,3),0), V(F(elId,3),1), V(F(elId,3),2));
                v4 = QVector3D(V(F(elId,4),0), V(F(elId,4),1), V(F(elId,4),2));
                v5 = QVector3D(V(F(elId,5),0), V(F(elId,5),1), V(F(elId,5),2));
                v6 = QVector3D(V(F(elId,6),0), V(F(elId,6),1), V(F(elId,6),2));
                v7 = QVector3D(V(F(elId,7),0), V(F(elId,7),1), V(F(elId,7),2));
                
                //Face 1
                //0,1,2
                n012 = QVector3D::normal(v0, v1, v2);
                
                //0,2,3
                n023 = QVector3D::normal(v0, v2, v3);
                
                //Face 2
                //7,6,5
                n765 = QVector3D::normal(v7, v6, v5);
                //7,5,4
                n754 = QVector3D::normal(v7, v5, v4);
                
                //Face 3
                //0,3,7
                n037 = QVector3D::normal(v0, v3, v7);
                //0,7,4
                n074 = QVector3D::normal(v0, v7, v4);
                
                //Face 4
                //6,2,1
                n621 = QVector3D::normal(v6, v2, v1);
                //5,6,1
                n561 = QVector3D::normal(v5, v6, v1);
                
                //Face 5
                //3,2,6
                n326 = QVector3D::normal(v3, v2, v6);
                //3,6,7
                n367 = QVector3D::normal(v3, v6, v7);
                
                //Face 6
                //1,0,4
                n104 = QVector3D::normal(v1, v0, v4);
                //1,4,5
                n145 = QVector3D::normal(v1, v4, v5);
                
                // Vector Normals
                n[0] = QVector3D(n023+n012+n037+n074+n104).normalized();
                n[1] = QVector3D(n012+n621+n561+n104+n145).normalized();
                n[2] = QVector3D(n012+n023+n621+n326).normalized();
                n[3] = QVector3D(n023+n037+n326+n367).normalized();
                n[4] = QVector3D(n754+n074+n104+n145).normalized();
                n[5] = QVector3D(n765+n754+n561+n145).normalized();
                n[6] = QVector3D(n765+n621+n561+n326+n367).normalized();
                n[7] = QVector3D(n765+n754+n037+n074+n367).normalized();
                
                for(unsigned int ii=0; ii< F.cols(); ++ii) {
                    rawVertexArray[idx++] = n[ii].x();
                    rawVertexArray[idx++] = n[ii].y();
                    rawVertexArray[idx++] = n[ii].z();
                }
                
                elId++;
                
                /*rawVertexArray[totalScalarV+3*F(elId,0)+0] = n0.x();
                 rawVertexArray[totalScalarV+3*F(elId,0)+1] = n0.y();
                 rawVertexArray[totalScalarV+3*F(elId,0)+2] = n0.z();
                 
                 rawVertexArray[totalScalarV+3*F(elId,1)+0] = n1.x();
                 rawVertexArray[totalScalarV+3*F(elId,1)+1] = n1.y();
                 rawVertexArray[totalScalarV+3*F(elId,1)+2] = n1.z();
                 
                 rawVertexArray[totalScalarV+3*F(elId,2)+0] = n2.x();
                 rawVertexArray[totalScalarV+3*F(elId,2)+1] = n2.y();
                 rawVertexArray[totalScalarV+3*F(elId,2)+2] = n2.z();
                 
                 rawVertexArray[totalScalarV+3*F(elId,3)+0] = n3.x();
                 rawVertexArray[totalScalarV+3*F(elId,3)+1] = n3.y();
                 rawVertexArray[totalScalarV+3*F(elId,3)+2] = n3.z();
                 elId++;*/
                
            }
            
            //idx += totalScalarV;
            
            //face colors
            elId = 0;
            
            //scale materials and use them to scale the colors
            //Let's just embrace lambdas
            auto ymMax = (*std::max_element(m_fem->getImpl().getElements().begin(), m_fem->getImpl().getElements().end(), [](auto a, auto b){ return a->getE() < b->getE(); }))->getE();
            auto ymMin = (*std::min_element(m_fem->getImpl().getElements().begin(),
                                            m_fem->getImpl().getElements().end(), [](auto a, auto b){ return a->getE() < b->getE(); }))->getE();
            
            
            std::cout<<"MAX/MIN "<<ymMax<<"/"<<ymMin<<"\n";
            double dym = (ymMax - ymMin);
            double color = 0.0;
            dym = (dym > 0 ? dym : 1.0);
            
            for(; elId < numElements; ++elId) {
                
                color =(m_fem->getImpl().getElement(elId)->getE() - ymMin)/dym;
                
                //std::cout<<"Color: "<<color<<"\n";
                
                for(unsigned int ii=0; ii<F.cols(); ++ii) {
                    rawVertexArray[idx++] = color*red.x();
                    rawVertexArray[idx++] = red.y();
                    rawVertexArray[idx++] = color*red.z();
                }
                
            }
            
            
            
            // Indices (12)
            QByteArray indexBufferData;
            indexBufferData.resize(36*numElements*sizeof(ushort)); //ah indices, for sweet sweet index rendering
            ushort *rawIndexArray = reinterpret_cast<ushort *>(indexBufferData.data());
            
            idx = 0;
            elId = 0;
            
            for(unsigned int idy = 0;idx < 36*numElements; ) {
                
                //Face 1
                //0,1,2
                rawIndexArray[idx++] = idy+0;
                rawIndexArray[idx++] = idy+1;
                rawIndexArray[idx++] = idy+2;
                
                //0,2,3
                rawIndexArray[idx++] = idy+0;
                rawIndexArray[idx++] = idy+2;
                rawIndexArray[idx++] = idy+3;
                
                
                //Face 2
                //7,6,5
                rawIndexArray[idx++] = idy+7;
                rawIndexArray[idx++] = idy+6;
                rawIndexArray[idx++] = idy+5;
                
                //7,5,4
                rawIndexArray[idx++] = idy+7;
                rawIndexArray[idx++] = idy+5;
                rawIndexArray[idx++] = idy+4;
                
                
                //Face 3
                //0,3,7
                rawIndexArray[idx++] = idy+0;
                rawIndexArray[idx++] = idy+3;
                rawIndexArray[idx++] = idy+7;
                
                //0,7,4
                rawIndexArray[idx++] = idy+0;
                rawIndexArray[idx++] = idy+7;
                rawIndexArray[idx++] = idy+4;
                
                
                //Face 4
                //6,2,1
                rawIndexArray[idx++] = idy+6;
                rawIndexArray[idx++] = idy+2;
                rawIndexArray[idx++] = idy+1;
                
                //5,6,1
                rawIndexArray[idx++] = idy+5;
                rawIndexArray[idx++] = idy+6;
                rawIndexArray[idx++] = idy+1;
                
                
                //Face 5
                //3,2,6
                rawIndexArray[idx++] = idy+3;
                rawIndexArray[idx++] = idy+2;
                rawIndexArray[idx++] = idy+6;
                
                //3,6,7
                rawIndexArray[idx++] = idy+3;
                rawIndexArray[idx++] = idy+6;
                rawIndexArray[idx++] = idy+7;
                
                
                //Face 6
                //1,0,4
                rawIndexArray[idx++] = idy+1;
                rawIndexArray[idx++] = idy+0;
                rawIndexArray[idx++] = idy+4;
                
                //1,4,5
                rawIndexArray[idx++] = idy+1;
                rawIndexArray[idx++] = idy+4;
                rawIndexArray[idx++] = idy+5;
            
                idy += 8;
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
            
             //6 faces of 6 points
            m_meshRenderer->setVertexCount(36*numElements);
            
            //transform
            m_transform = new Qt3DCore::QTransform();
            m_transform->setScale3D(QVector3D(1,1,1));
            m_transform->setRotationX(0.0);
            m_transform->setRotationY(0.0);
            m_transform->setRotationZ(0.0);
            m_transform->setTranslation(QVector3D(0,0,0));
            
            //Qt3DExtras::QPhongMaterial *material = new Qt3DExtras::QPhongMaterial();
            //Qt3DRender::QMaterial *material = new Qt3DExtras::QPerVertexColorMaterial();
            
            //material->setDiffuse(QColor(220.0, 0.0,0.0));
            
            m_entity = new Qt3DCore::QEntity(root);
            m_entity->addComponent(m_meshRenderer);
            m_entity->addComponent(m_transform);
            m_entity->addComponent(wireframeMaterial);
            
            //setup update array
            
            return m_entity;
            
        }
        
        void update(State<DataType> &state) const {
            
            //reset mesh positions and add in the state (which I think is already setup in the right order ?)
            
            
            //1.) Grab raw data
            float *rawData = const_cast<float *>(reinterpret_cast<const float *>(m_updateData.data()));
            float *rawOriginalData = const_cast<float *>(reinterpret_cast<const float *>(m_vertexBufferData.data()));
            
            //2.) Update using mesh data + state which for linear FEM is the displacement of the node
            Eigen::Map<Eigen::VectorXd> pos = mapDOFEigen(m_fem->getQ(), state);
            
            unsigned int vertexId, idx;
            idx = 0;
            
            for(unsigned int elId=0;  elId < m_fem->getImpl().getF().rows(); ++elId) {
                
                for(unsigned int ii=0; ii<m_fem->getImpl().getF().cols(); ++ii) {
                    vertexId = m_fem->getImpl().getF()(elId, ii);
                    rawData[idx] = rawOriginalData[idx] + pos(m_fem->getImpl().getQ()[vertexId].getLocalId());
                    idx++;
                    rawData[idx] = rawOriginalData[idx] + pos(m_fem->getImpl().getQ()[vertexId].getLocalId()+1);
                    idx++;
                    rawData[idx] = rawOriginalData[idx] + pos(m_fem->getImpl().getQ()[vertexId].getLocalId()+2);
                    idx++;
                }
                
                /*vertexId = m_fem->getImpl().getF()(ii, 0);
                rawData[idx] = rawOriginalData[idx] + pos(m_fem->getImpl().getQ()[vertexId].getGlobalId());
                idx++;
                rawData[idx] = rawOriginalData[idx] + pos(m_fem->getImpl().getQ()[vertexId].getGlobalId()+1);
                idx++;
                rawData[idx] = rawOriginalData[idx] + pos(m_fem->getImpl().getQ()[vertexId].getGlobalId()+2);
                idx++;
                
                //vertex 2
                vertexId = m_fem->getImpl().getF()(ii, 1);
                rawData[idx] = rawOriginalData[idx] + pos(m_fem->getImpl().getQ()[vertexId].getGlobalId());
                idx++;
                rawData[idx] = rawOriginalData[idx] + pos(m_fem->getImpl().getQ()[vertexId].getGlobalId()+1);
                idx++;
                rawData[idx] = rawOriginalData[idx] + pos(m_fem->getImpl().getQ()[vertexId].getGlobalId()+2);
                idx++;
                
                //vertex 3
                vertexId = m_fem->getImpl().getF()(ii, 2);
                rawData[idx] = rawOriginalData[idx] + pos(m_fem->getImpl().getQ()[vertexId].getGlobalId());
                idx++;
                rawData[idx] = rawOriginalData[idx] + pos(m_fem->getImpl().getQ()[vertexId].getGlobalId()+1);
                idx++;
                rawData[idx] = rawOriginalData[idx] + pos(m_fem->getImpl().getQ()[vertexId].getGlobalId()+2);
                idx++;
                
                //vertex 4
                vertexId = m_fem->getImpl().getF()(ii, 3);
                rawData[idx] = rawOriginalData[idx] + pos(m_fem->getImpl().getQ()[vertexId].getGlobalId());
                idx++;
                rawData[idx] = rawOriginalData[idx] + pos(m_fem->getImpl().getQ()[vertexId].getGlobalId()+1);
                idx++;
                rawData[idx] = rawOriginalData[idx] + pos(m_fem->getImpl().getQ()[vertexId].getGlobalId()+2);
                idx++;*/
                
            }
            
            //Update
            m_vertexDataBuffer->updateData(0, m_updateData);
            
        }
        
    protected:
        Qt3DCore::QEntity *m_entity;
        Qt3DRender::QGeometryRenderer *m_meshRenderer;
        Qt3DRender::QGeometry *m_meshGeometry;
        FEMSystem *m_fem;
        Qt3DCore::QTransform *m_transform;
        QByteArray m_updateData;
        QByteArray m_vertexBufferData;
        Qt3DRender::QBuffer *m_vertexDataBuffer;
        
    private:
    };

}

#endif /* RenderableFEMHex_h */
