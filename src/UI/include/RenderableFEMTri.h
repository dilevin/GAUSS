//
//  RenderableFEMTri.h
//  Gauss
//
//  Created by David Levin on 1/5/18.
//
//

#ifndef RenderableFEMTri_h
#define RenderableFEMTri_h

template<typename DataType, template <typename D, typename E> class QuadratureRuleT, template <typename H, typename I> class QuadratureRuleU, template <typename F, typename G> class EnergyPotential>
class Renderable<PhysicalSystem<DataType, FEM::PhysicalSystemFEMImpl<DataType,
FEM::ElementBase<DataType, 3, QuadratureRuleT, QuadratureRuleU,
FEM::EnergyKineticNonLumped,
EnergyPotential,
FEM::BodyForceGravity,
FEM::ShapeFunctionPlaneLinear> > > >: public Renderable<DataType>
{
public:
    
    using FEMSystem = PhysicalSystem<DataType, FEM::PhysicalSystemFEMImpl<DataType,
    FEM::ElementBase<DataType, 3, QuadratureRuleT, QuadratureRuleU,
    FEM::EnergyKineticNonLumped,
    EnergyPotential,
    FEM::BodyForceGravity,
    FEM::ShapeFunctionPlaneLinear> > >;
    
    
    Renderable(FEMSystem *fem) : Renderable<DataType>() {
        std::cout<<"Renderable Tet FEM \n";
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
        m_tetRenderer = new Qt3DRender::QGeometryRenderer;
        
        //custom mesh renderer takes in geometry as well as information about offset per primitive and type of primitive to render
        m_tetGeometry = new Qt3DRender::QGeometry(m_tetRenderer);
        
        m_vertexDataBuffer = new Qt3DRender::QBuffer(Qt3DRender::QBuffer::VertexBuffer, m_tetGeometry);
        
        Qt3DRender::QBuffer *indexDataBuffer = new Qt3DRender::QBuffer(Qt3DRender::QBuffer::IndexBuffer, m_tetGeometry);
        
        
        //long totalScalarV = m_fem->getImpl().getV().rows()*m_fem->getImpl().getV().cols();
        
        //total numebr of verts = numVetsPerElement*numElements
        long totalScalarV = 3*m_fem->getImpl().getF().cols()*m_fem->getImpl().getF().rows();
        long numTets = m_fem->getImpl().getF().rows();
        
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
        int tetId = 0;
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
        
        for(; idx<totalScalarV; ++tetId) {
            
            //vertex 1
            vertexId = F(tetId, 0);
            rawUpdatePositionArray[idx] = m_fem->getImpl().getV()(vertexId,0);
            rawVertexArray[idx++] = m_fem->getImpl().getV()(vertexId,0);
            rawUpdatePositionArray[idx] = m_fem->getImpl().getV()(vertexId,1);
            rawVertexArray[idx++] = m_fem->getImpl().getV()(vertexId,1);
            rawUpdatePositionArray[idx] = 0;
            rawVertexArray[idx++] = 0;
            
            
            //vertex 2
            vertexId = F(tetId, 1);
            rawUpdatePositionArray[idx] = m_fem->getImpl().getV()(vertexId,0);
            rawVertexArray[idx++] = m_fem->getImpl().getV()(vertexId,0);
            rawUpdatePositionArray[idx] = m_fem->getImpl().getV()(vertexId,1);
            rawVertexArray[idx++] = m_fem->getImpl().getV()(vertexId,1);
            rawUpdatePositionArray[idx] = 0;
            rawVertexArray[idx++] = 0;
            
            //vertex 3
            vertexId = F(tetId, 2);
            rawUpdatePositionArray[idx] = m_fem->getImpl().getV()(vertexId,0);
            rawVertexArray[idx++] = m_fem->getImpl().getV()(vertexId,0);
            rawUpdatePositionArray[idx] = m_fem->getImpl().getV()(vertexId,1);
            rawVertexArray[idx++] = m_fem->getImpl().getV()(vertexId,1);
            rawUpdatePositionArray[idx] = 0;
            rawVertexArray[idx++] = 0;
            
            
        }
        
        QVector3D v0, v1, v2, v3;
        
        
        //really lazy vertex normals for now
        tetId = 0;
        
        for(; tetId < numTets; ) {
            
         
            // Vector Normals
            
            rawVertexArray[idx++] = 0;
            rawVertexArray[idx++] = 0;
            rawVertexArray[idx++] = 1;
            
            rawVertexArray[idx++] = 0;
            rawVertexArray[idx++] = 0;
            rawVertexArray[idx++] = 1;
            
            rawVertexArray[idx++] = 0;
            rawVertexArray[idx++] = 0;
            rawVertexArray[idx++] = 1;
            
            tetId++;
            
        }
        
        //idx += totalScalarV;
        
        //face colors
        tetId = 0;
        
        //scale materials and use them to scale the colors
        //Let's just embrace lambdas
        auto ymMax = (*std::max_element(m_fem->getImpl().getElements().begin(),
                                        m_fem->getImpl().getElements().begin(), [](auto a, auto b){ return a->getE() < b->getE(); }))->getE();
        auto ymMin = (*std::min_element(m_fem->getImpl().getElements().begin(),
                                        m_fem->getImpl().getElements().begin(), [](auto a, auto b){ return a->getE() < b->getE(); }))->getE();
        
        
        double dym = (ymMax - ymMin);
        double color = 0.0;
        dym = (dym > 0 ? dym : 1.0);
        
        for(; tetId < numTets; ++tetId) {
            
            color =(m_fem->getImpl().getElement(tetId)->getE() - ymMin)/dym;
            color = 0.0;
            rawVertexArray[idx++] = color*red.x();
            rawVertexArray[idx++] = red.y();
            rawVertexArray[idx++] = color*red.z();
            
            rawVertexArray[idx++] = color*red.x();
            rawVertexArray[idx++] = red.y();
            rawVertexArray[idx++] = color*red.z();
            
            rawVertexArray[idx++] = color*red.x();
            rawVertexArray[idx++] = red.y();
            rawVertexArray[idx++] = color*red.z();
            
        }
        
        
        // Indices (12)
        QByteArray indexBufferData;
        indexBufferData.resize(12*numTets*sizeof(ushort)); //ah indices, for sweet sweet index rendering
        ushort *rawIndexArray = reinterpret_cast<ushort *>(indexBufferData.data());
        
        idx = 0;
        tetId = 0;
        for(unsigned int idy = 0;idx < 3*numTets; ) {
            // Front
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
        
        m_tetGeometry->addAttribute(positionAttribute);
        m_tetGeometry->addAttribute(normalAttribute);
        m_tetGeometry->addAttribute(colorAttribute);
        m_tetGeometry->addAttribute(indexAttribute);
        
        m_tetRenderer->setInstanceCount(1);
        m_tetRenderer->setIndexOffset(0);
        m_tetRenderer->setFirstInstance(0);
        m_tetRenderer->setPrimitiveType(Qt3DRender::QGeometryRenderer::Triangles);
        m_tetRenderer->setGeometry(m_tetGeometry);
        
        // 4 faces of 3 points
        m_tetRenderer->setVertexCount(3*numTets);
        
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
        m_entity->addComponent(m_tetRenderer);
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
        for(unsigned int ii=0;  ii < m_fem->getImpl().getF().rows(); ++ii) {
            
            vertexId = m_fem->getImpl().getF()(ii, 0);
            rawData[idx] = rawOriginalData[idx];// + pos(m_fem->getImpl().getQ()[vertexId].getGlobalId());
            idx++;
            rawData[idx] = rawOriginalData[idx] + pos(m_fem->getImpl().getQ()[vertexId].getGlobalId()+1);
            idx++;
            rawData[idx] = 0;
            idx++;
            
            //vertex 2
            vertexId = m_fem->getImpl().getF()(ii, 1);
            rawData[idx] = rawOriginalData[idx] + pos(m_fem->getImpl().getQ()[vertexId].getGlobalId());
            idx++;
            rawData[idx] = rawOriginalData[idx] + pos(m_fem->getImpl().getQ()[vertexId].getGlobalId()+1);
            idx++;
            rawData[idx] = 0;
            idx++;
            
            //vertex 3
            vertexId = m_fem->getImpl().getF()(ii, 2);
            rawData[idx] = rawOriginalData[idx] + pos(m_fem->getImpl().getQ()[vertexId].getGlobalId());
            idx++;
            rawData[idx] = rawOriginalData[idx] + pos(m_fem->getImpl().getQ()[vertexId].getGlobalId()+1);
            idx++;
            rawData[idx] = 0;
            idx++;
            
        }
    
    
        //Update
        m_vertexDataBuffer->updateData(0, m_updateData);
        
    }
    
protected:
    Qt3DCore::QEntity *m_entity;
    Qt3DRender::QGeometryRenderer *m_tetRenderer;
    Qt3DRender::QGeometry *m_tetGeometry;
    FEMSystem *m_fem;
    Qt3DCore::QTransform *m_transform;
    QByteArray m_updateData;
    QByteArray m_vertexBufferData;
    Qt3DRender::QBuffer *m_vertexDataBuffer;
    
private:
};


#endif /* RenderableFEMTri_h */
