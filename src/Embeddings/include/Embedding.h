#ifndef _EMBEDDINGS_H
#define _EMBEDDINGS_H
#include <GaussIncludes.h>
#include <FEMIncludes.h>

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
                getShapeFunctionMatrix(N,V, fem);
                m_N = (*N);
                
            }
            
            inline decltype(auto) getGeometry() {
                
                return std::make_pair(std::ref(m_V), std::ref(m_F));
            }
            
            inline const auto getPosition(const PhysicalSystemImpl &fem, const State<DataType> &state, unsigned int vertexId) const {
                
                return m_V.row(vertexId).transpose() + m_N.block(3*vertexId, 0, 3, m_N.cols())*mapDOFEigen(fem.getQ(), state);
            }
            
                                               
            inline void getDPDQ(const PhysicalSystemImpl &fem, const State<DataType> &state, unsigned int vertexId) {
                std::cout<<"DPDQ Not implemented yet for embedded mesh \n";
                exit(1);
            }
            
            inline void getVelocity(const PhysicalSystemImpl &fem, const State<DataType> &state, unsigned int vertexId) {
                return m_N.block(3*vertexId, 0, 3, m_N.cols())*mapDOFEigen(fem.getQDot(), state);
            }
            
            inline void getDVDQ(const PhysicalSystemImpl &fem, const State<DataType> &state, unsigned int vertexId) {
                std::cout<<"DVDQ Not implemented yet for embedded mesh \n";
                exit(1);
            }
            
            
        protected:
            
           Eigen::MatrixXx<DataType> m_V;
           Eigen::MatrixXi  m_F;
           Eigen::SparseMatrix<DataType, Eigen::RowMajor> m_N; //embedding matrix
            
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
            
            //contructor
            template<typename ...Params>
            PhysicalSystemEmbeddedMeshImpl(Eigen::MatrixXx<DataType> &V, Eigen::MatrixXi &F, Params ...params) : Embedding::PhysicalSystemImpl(params...), m_embedding(V,F, *this) { }
            
            //Geometry stuff
            template<typename ...Params>
            inline decltype(auto) getPosition(const State<DataType> &state, Params &...params) const {
                return m_embedding.getPosition(*this, state, params...);
            }
            
            template<typename Vector, typename ...Params>
            inline auto getDPDQ(Vector &x, Params &... params) {
                return m_embedding.getDPDQ(*this, x, params...);
            }
            
            template<typename  Vector, typename ...Params>
            inline auto getVelocity(Vector &x, Params &...params) {
                return m_embedding.getVelocity(*this, x, params...);
            }
            
            template<typename  Vector, typename ...Params>
            inline auto getDVDQ(Vector &x, Params &...params) {
                return m_embedding.getDVDQ(*this, x, params...);
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
