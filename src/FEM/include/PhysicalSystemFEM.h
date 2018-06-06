//Design Notes:
//Physical System:
//MultiVector<ElementTypes>
//MultiVector<DOFType>
//Store DOFs in continguous memory
//Each element
//  Quadrature
//      Energy(ShapeFunction, Energy, position) -> Energy(ShapeFunction(position));
//  ShapeFunction + Kinematics (DOF + gradients I might need ....)
//      Gradient
//      Hessian (maybe, I'm not sure)
//      How to represent DOFS ? Pointers to list from Physical System ? Has to be, since elements like forces are   connections
//  Energy
//      //Gradient
//      //Hessian
//How should lists of DOFS work ?
//  Container<DOFs>
#ifndef PHYSICALSYSTEMFEM_H
#define PHYSICALSYSTEMFEM_H

#include <vector>
#include <DOFParticle.h>
#include <DOFList.h>
#include <UtilitiesEigen.h>

namespace Gauss {
    namespace FEM {
        template<typename DataType, typename ElementType>
        class PhysicalSystemFEMImpl
        {
            
        public:
            //temporary global indices until I update the state to give these to me
            //automatically
            PhysicalSystemFEMImpl(const Eigen::Ref<Eigen::MatrixXd> &V, const Eigen::Ref<Eigen::MatrixXi> &F) : m_q(V.rows()), m_qDot(V.rows()) {
                
                m_V = V;
                m_F = F;
                m_numVerts = m_V.rows();
                m_numElements = m_F.rows();
                
                
                assert(m_V.cols() == 3); //3D only for now
                
                //initialize all the elements
                Eigen::MatrixXi Felement;
                std::array<DOFBase<DataType,0> *, ElementType::numDOFs()> qDOFArray;
                std::array<DOFBase<DataType,1> *, ElementType::numDOFs()> qDotDOFArray;
                for(unsigned int iel=0; iel < m_numElements; iel++) {
                    
                    for(unsigned int idof=0;idof < ElementType::numDOFs(); ++idof) {
                        qDOFArray[idof] = &m_q[F(iel,idof)];
                        qDotDOFArray[idof] = &m_qDot[F(iel,idof)];
                    }
                    
                    Felement = m_F.row(iel);
                   
                    m_elements.push_back(
                       new ElementType(m_V,Felement, qDOFArray, qDotDOFArray)
                    );
                }
                
                
            }
            
            ~PhysicalSystemFEMImpl() {
                
            }
            
            DataType getEnergy(const State<DataType> &state) const {
                
                double energy = 0.0;
                for(auto &element : m_elements) {
                    energy += element->getEnergy(state);
                }
                
                return energy;
            }

            DataType getBodyForceEnergy(const State<DataType> &state) const {
                DataType energy = 0.0;
                for(auto &element : m_elements) {
                    energy += element->getBodyForceWork(state);
                }
                
                return energy;
            }
            
            DataType getStrainEnergy(const State<DataType> &state) const {
                
                DataType energy = 0.0;
                for(auto &element : m_elements) {
                    energy += element->getStrainEnergy(state);
                }
                
                return energy;
            }

            decltype(auto) getStrainEnergyPerElement(const State<DataType> &state) const {
                Eigen::VectorXx<DataType> energyPerElement(m_elements.size());

                for(int i=0; i < m_elements.size(); i++) {
                    energyPerElement[i] = m_elements[i]->getStrainEnergy(state);                
                }
                
                return energyPerElement;
            }
            
            template<typename Assembler>
            inline void getMassMatrix(Assembler &assembler, const State<DataType> &state) const {
                //call the assembler on all elements
                forLoop<IsParallel<Assembler>::value>(m_elements, assembler, [&](auto &assemble, auto &element) {
                    element->getMassMatrix(assemble,state);
                });
            }
            
            template<typename Assembler>
            inline void getStiffnessMatrix(Assembler &assembler, const State<DataType> &state) const {
                
                forLoop<IsParallel<Assembler>::value>(m_elements, assembler, [&](auto &assemble, auto &element) {
                    element->getStiffnessMatrix(assemble, state);
                });
            }
            
            template<typename Assembler>
            inline void getForce(Assembler &assembler, const State<DataType> &state) const {
        
                forLoop<IsParallel<Assembler>::value>(m_elements, assembler, [&](auto &assemble, auto &element) {
                    element->getForce(assemble, state);
                });
            }
            
            template<typename Assembler>
            inline void getInternalForce(Assembler &assembler, const State<DataType> &state) const {
  
                forLoop<IsParallel<Assembler>::value>(m_elements, assembler, [&](auto &assemble, auto &element) {
                    element->getInternalForce(assemble, state);
                });
            }
            
            template<typename Assembler>
            inline void getBodyForce(Assembler &assembler, const State<DataType> &state) const {
                forLoop<IsParallel<Assembler>::value>(m_elements, assembler, [&](auto &assemble, auto &element) {
                    element->getBodyForce(assemble, state);
                });
            }
            
            inline unsigned int getNumElements() { return m_elements.size(); }
            
            inline ElementType * getElement(unsigned int i) {
                assert(i < m_elements.size());
                return m_elements[i];
                
            }
            
            inline std::vector<ElementType *> & getElements() { return m_elements; }
            inline const std::vector<ElementType *> & getElements() const { return m_elements; }
            
            inline const ElementType * getElement(unsigned int i) const {
                assert(i < m_elements.size());
                return m_elements[i];
                
            }
            
            inline auto & getQ() { return m_q; }
            inline const auto & getQ() const { return m_q; }
            
            inline auto & getQDot() { return m_qDot; }
            inline const auto & getQDot() const { return m_qDot; }

            //get function supporting a vertex (these return arrays in order to slot directly into assemblers)
            inline decltype(auto) getQ(unsigned int vertexId) const {
                std::array<const DOFBase<DataType,0> *,1> toReturn = {{&m_q[vertexId]}};
                return toReturn;
            }
           
            inline decltype(auto) getQDot(unsigned int vertexId) const {
                std::array<const DOFBase<DataType,1> *,1> toReturn = {{&m_qDot[vertexId]}};
                return toReturn;
            }
            
            
            template<typename Vector>
            inline decltype(auto) getQ(Vector &x, unsigned int elementId) const {
                std::cout<<"Error not implemented \n";
                exit(0);
                std::array<const DOFBase<DataType,0> *, 1> toReturn = {{&m_q[elementId]}};
                return toReturn;
            }
            
            template<typename Vector>
            inline decltype(auto) getQDot(Vector &x, unsigned int elementId) const {
                std::cout<<"Error not implemented \n";
                exit(0);
                std::array<const DOFBase<DataType,1> *,1> toReturn = {{&m_qDot[elementId]}};
                return toReturn;
            }
            
            inline auto & getV() { return m_V; }
            inline auto & getF() { return m_F; }
            
            inline const auto & getV() const { return m_V; }
            inline const auto & getF() const { return m_F; }
            
           
            //methods for getting current positions and position Jacobians for this system
            //Per-Vertex
            inline const auto getPosition(const State<DataType> &state, unsigned int vertexId) const {
                return getV().row(vertexId).transpose() + mapDOFEigen(m_q[vertexId], state);
            }
            
            inline const auto getVelocity(const State<DataType> &state, unsigned int vertexId) const {
                return mapDOFEigen(m_qDot[vertexId], state);
            }
            
            inline const auto getDPDQ(const State<DataType> &state, unsigned int vertexId) const {
                return Eigen::Matrix33x<DataType>::Identity();
            }
            
            inline const auto getDPDQ(const State<DataType> &state, unsigned int elementId, const Eigen::Vector3x<DataType> &pos) const {
                exit(0);
                return Eigen::Matrix33x<DataType>::Identity();
            }
            
            //want these for elements as well (i.e take an element indec and a point in space and return the right value)
            
            
            inline auto getGeometry() { return std::make_pair(std::ref(m_V), std::ref(m_F)); }
            
        protected:
            
            //Mesh
            Eigen::MatrixXd m_V;
            Eigen::MatrixXi m_F;
            long m_numVerts;
            long m_numElements;
           
            DOFList<DataType, DOFParticle, 0> m_q;
            DOFList<DataType, DOFParticle, 1> m_qDot;
            std::vector<ElementType *> m_elements;
            //DataType m_mass; //mass of particle
            //DOFParticle<DataType,0> m_x;
            //DOFParticle<DataType,1> m_xDot;
            
        private:
  
        };
 
        template<typename DataType, template <typename A>  class ElementType>
        using  PhysicalSystemFEM = PhysicalSystem<DataType, PhysicalSystemFEMImpl<DataType, ElementType<DataType> > >;

    }
}

#endif
