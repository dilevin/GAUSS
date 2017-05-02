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
            PhysicalSystemFEMImpl(Eigen::MatrixXd &V, Eigen::MatrixXi &F) : m_q(V.rows()), m_qDot(V.rows()) {
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
                    qDOFArray[0] = &m_q[F(iel,0)];
                    qDOFArray[1] = &m_q[F(iel,1)];
                    qDOFArray[2] = &m_q[F(iel,2)];
                    qDOFArray[3] = &m_q[F(iel,3)];
                    qDotDOFArray[0] = &m_qDot[F(iel,0)];
                    qDotDOFArray[1] = &m_qDot[F(iel,1)];
                    qDotDOFArray[2] = &m_qDot[F(iel,2)];
                    qDotDOFArray[3] = &m_qDot[F(iel,3)];
                    
                    Felement = m_F.row(iel);
                   
                    m_elements.push_back(
                       new ElementType(m_V,Felement, qDOFArray, qDotDOFArray)
                    );
                }
                
            }
            
            ~PhysicalSystemFEMImpl() {
                
            }
            
            template<typename Assembler>
            inline void getMassMatrix(Assembler &assembler, const State<DataType> &state) const {
                //call the assembler on all elements
                for(auto element : m_elements) {
                    element->getMassMatrix(assembler, state);
                }
            }
            
            template<typename Assembler>
            inline void getStiffnessMatrix(Assembler &assembler, const State<DataType> &state) const {
                
                for(auto element : m_elements) {
                    element->getStiffnessMatrix(assembler, state);
                }
            }
            
            template<typename Assembler>
            inline void getForce(Assembler &assembler, const State<DataType> &state) const {
                for(auto element : m_elements) {
                    element->getForce(assembler, state);
                }
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
            
            auto & getQ() { return m_q; }
            const auto & getQ() const { return m_q; }
            
            auto & getQDot() { return m_qDot; }
            const auto & getQDot() const { return m_qDot; }
            
            inline auto & getV() { return m_V; }
            inline auto & getF() { return m_F; }
            
            inline const auto & getV() const { return m_V; }
            inline const auto & getF() const { return m_F; }
            
            
        protected:
            
            //Mesh
            Eigen::MatrixXd m_V;
            Eigen::MatrixXi m_F;
            long m_numVerts;
            long m_numElements;
            
            DOFList<DataType, ParticleSystem::DOFParticle, 0> m_q;
            DOFList<DataType, ParticleSystem::DOFParticle, 1> m_qDot;
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
