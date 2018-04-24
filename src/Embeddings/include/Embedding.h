#ifndef _EMBEDDINGS_H
#define _EMBEDDINGS_H
#include <GaussIncludes.h>
#include <FEMIncludes.h>
#include <igl/barycentric_coordinates.h>

namespace Gauss {
    namespace Embeddings {
        
        template<typename DataType, typename EmbeddingSystem>
        class EmbeddingFunction {
          
        public:
            using PhysicalSystemImpl = EmbeddingSystem;
            
            template<typename ...Params>
            EmbeddingFunction(Params ...params) {
               std::cout<<"Generic Embedding \n";
            }
            
            inline void getGeometry() {  }
            
            inline void getPosition() {  }
            
            inline void getDPDQ() { }
            
            inline void getVelocity() { }
            
            inline void getDVDQ() { }
            
            
        protected:
            
           
        private:
        
        };
        
        template<typename DataType,
                unsigned int NUM,
                template <typename D, typename E> class QuadratureRuleT,
                template <typename H, typename I> class QuadratureRuleU,
                template <typename J, typename K> class EnergyKinetic,
                template <typename L, typename M> class EnergyPotential,
                template <typename N, typename O> class BodyForce,
                template <typename P> class ShapeFunction
                >
        class EmbeddingFunction<DataType,  FEM::PhysicalSystemFEMImpl<DataType,
                                           FEM::ElementBase<DataType, NUM,
                                           QuadratureRuleT,
                                           QuadratureRuleU,
                                           EnergyKinetic,
                                           EnergyPotential,
                                           BodyForce,
                                           ShapeFunction> > >{
            
        public:
              using PhysicalSystemImpl =  FEM::PhysicalSystemFEMImpl<DataType,
                                               FEM::ElementBase<DataType, NUM,
                                               QuadratureRuleT,
                                               QuadratureRuleU,
                                               EnergyKinetic,
                                               EnergyPotential,
                                               BodyForce,
                                               ShapeFunction> > ;
            
            template<typename ...Params>
            EmbeddingFunction(Eigen::MatrixXx<DataType> &V, Eigen::MatrixXi &F, PhysicalSystemImpl &fem) {
                std::cout<<"FEM Embedding \n";
                
                m_V = V;
                m_F = F;
                
                AssemblerEigenSparseMatrix<DataType> N;
                //element[i] is a n-vector that stores the index of the element containing the ith vertex in the embedded mesh
                getShapeFunctionMatrix(N,m_elements,V, fem);
            
                Eigen::Vector3x<DataType> vertex = m_V.row(0);
                
                unsigned int numCols = fem.getElements()[0]->N(vertex.data()).cols();
                unsigned int el;
                
                m_N.resize(3*m_V.rows(), numCols);
                for(unsigned int ii=0; ii<m_elements.size(); ++ii) {
                    el = m_elements[ii];
                    vertex = m_V.row(ii);
                    
                    auto Jmat = fem.getElements()[el]->N(vertex.data());
                    
                    m_N.block(3*ii, 0, 3, numCols) = Jmat;
                }
                
                
            }
            
            inline decltype(auto) getGeometry() {
                
                return std::make_pair(std::ref(m_V), std::ref(m_F));
            }
            

            //per vertex accessors
            inline Eigen::Vector3x<DataType> getPosition(const PhysicalSystemImpl &fem, const State<DataType> &state, unsigned int vertexId) const {
                //return m_V.row(vertexId).transpose() + m_N.block(3*vertexId, 0, 3, m_N.cols())*mapDOFEigen(fem.getQ(), state);
                return m_V.row(vertexId).transpose() + m_N.block(3*vertexId, 0, 3, m_N.cols())*fem.getElement(m_elements[vertexId])->q(state);
            }
            
            //this is also surface J
            inline Eigen::MatrixXx<DataType>  getDPDQ(const PhysicalSystemImpl &fem, const State<DataType> &state, unsigned int vertexId) const {
                return m_N.block(3*vertexId, 0, 3, m_N.cols());
            }
            
            inline Eigen::Vector3x<DataType> getVelocity(const PhysicalSystemImpl &fem, const State<DataType> &state, unsigned int vertexId) {
                return m_N.block(3*vertexId, 0, 3, m_N.cols())*fem.getElement(m_elements[vertexId])->qDot(state);
            }
            
            //per spatial point accessors (nothing implemented for these yet)
            inline Eigen::MatrixXx<DataType> getDPDQ(const PhysicalSystemImpl &fem, const State<DataType> &state, unsigned int elementId, const Eigen::Vector3x<DataType> &x) const {
               
                //Eigen::MatrixXx<DataType> J;
                //J.resize(3, 3*m_N.cols());
                Eigen::MatrixXx<DataType> J;
                J.resize(3, 3*m_N.cols());
                
                J.block(0,0, 3, m_N.cols()) = getDPDQ(fem, state, m_F(elementId, 0));
                J.block(0,m_N.cols(), 3, m_N.cols()) = getDPDQ(fem, state, m_F(elementId, 1));
                J.block(0,2*m_N.cols(), 3, m_N.cols()) = getDPDQ(fem, state, m_F(elementId, 2));
                
                return J;
                
           }
            
           template<typename Vector>
           inline decltype(auto) getQ(const PhysicalSystemImpl &fem, Vector &x, unsigned int elementId)  const {
               
               std::array<DOFBase<DataType> *, 3*NUM> q;
               unsigned int kk=0;
               for(unsigned int ii=0; ii< 3; ++ii) {
                   auto e = fem.getElement(m_elements[m_F(elementId, ii)])->q();
                   for(unsigned int jj=0 ;jj< e.size(); ++jj) {
                       q[kk] = e[jj];
                       kk++;
                   }
               }
               return q;
           }
                                               
           inline decltype(auto) getQ(const PhysicalSystemImpl &fem, unsigned int vertexId) const {

               return fem.getElement(m_elements[vertexId])->q();
           }
                                               
           inline  decltype(auto) getQDot(const PhysicalSystemImpl &fem, unsigned int vertexId) const  {
               
               return fem.getElement(m_elements[vertexId])->qDot();
           }
                
        protected:
            
           Eigen::MatrixXx<DataType> m_V;
           Eigen::MatrixXi  m_F;
           Eigen::MatrixXd m_N;
           //m_elements[i] is a n-vector that stores the index of the element containing the ith vertex in the embedded mesh
           Eigen::VectorXi m_elements;
           
                                               
        private:
            
        };
        
        template<typename DataType, typename System>
        class  PhysicalSystemEmbeddedMeshImpl : public EmbeddingFunction<DataType, typename System::ImplType>::PhysicalSystemImpl {
        
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
            
            //embedding dense matrix size
            
            //contructor
            template<typename ...Params>
            PhysicalSystemEmbeddedMeshImpl(Eigen::MatrixXx<DataType> &V, Eigen::MatrixXi &F, Params ...params) : Embedding::PhysicalSystemImpl(params...), m_embedding(V,F, *this) { }
            
            // First two parameters V and F are for the embedded high res mesh in the constructor for m_embedding, line 68 of Embedding.h
            // the next two parameters are V and F from the coarse mesh, packed in ...Params and unpacked in ...params, which is than
            // passed to the constructor for Embedding::PhysicalSystemImpl, which is actually Embedding::FEM::PhysicalSystemFEMImpl, and its constructor is
            // on line 35 of PhysicalSystemFEM.h
            
            //Geometry stuff
            template<typename ...Params>
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
            
            inline decltype(auto) getGeometry() { return m_embedding.getGeometry(); }
            
        protected:
            Embedding m_embedding;
            
            
        private:

            
        };
        
        template<typename DataType, typename System>
        using PhysicalSystemEmbeddedMesh = PhysicalSystem<DataType, PhysicalSystemEmbeddedMeshImpl<DataType, System> >;
    }
}
#endif
