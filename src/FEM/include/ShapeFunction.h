//
//  ShapeFunction.h
//  Gauss
//
//  Created by David Levin on 3/13/17.
//
//

#ifndef ShapeFunction_h
#define ShapeFunction_h

#include <array>
#include <Utilities.h>

//How this should actually work:
// Quadrature : Energy : Implementation: Shape Functions
//Shape functions just stores nodal, scalar functiosn and can compute their derivatives + derivative
//Kinematics: Handles getting deformation gradient and all that stuff
//Everything else same as I have it
//But due to laziness for now I'll concatanate Kinematics and ShapeFunctions into a single class
//Can split them up later when its convenient
//Let's use V,F for everythign so we can just hand these around to initialize things
namespace Gauss {
    namespace FEM {
        //Blank Shape Function
        template<typename DataType>
        class ShapeFunctionNone {
        public:
            template <typename Derived,typename DOFList>
            ShapeFunctionNone(Eigen::MatrixXd &V, Eigen::MatrixBase<Derived> F, DOFList &dofList) { std::cout<<"Do nothing Shape Function \n"; }
        };
        
        //Linear Tetrahedral Shape Function
        template<typename DataType>
        class ShapeFunctionLinearTet {
        public:
            
            template<typename QDOFList, typename QDotDOFList>
            ShapeFunctionLinearTet(Eigen::MatrixXd &V, Eigen::MatrixXi &F, QDOFList &qDOFList, QDotDOFList &qDotDOFList) : m_T(3,3) {
                //build up stuff I need for barycentric coordinates
              
                m_T << (Vert(1,0) - Vert(0,0)), (Vert(2,0) - Vert(0,0)), (Vert(3,0) - Vert(0,0)),
                       (Vert(1,1) - Vert(0,1)), (Vert(2,1) - Vert(0,1)), (Vert(3,1) - Vert(0,1)),
                       (Vert(1,2) - Vert(0,2)), (Vert(2,2) - Vert(0,2)), (Vert(3,2) - Vert(0,2));
                
                
                //invert
                
                m_x3 << Vert(0,0),
                        Vert(0,1),
                        Vert(0,2);
                
                m_T = m_T.inverse().eval();
                
                //for the time being assume things come as a stack (qdofs and qdotdofs)
                m_qDofs[0] = qDOFList[0];
                m_qDofs[1] = qDOFList[1];
                m_qDofs[2] = qDOFList[2];
                m_qDofs[3] = qDOFList[3];
                
                m_qDotDofs[0] = qDotDOFList[0];
                m_qDotDofs[1] = qDotDOFList[1];
                m_qDotDofs[2] = qDotDOFList[2];
                m_qDotDofs[3] = qDotDOFList[3];
            }
            
            //drop my gets for this just because my hands get tired of typing it all the time
            template<unsigned int Vertex>
            inline double phi(double *x) {
                assert(Vertex > 0);
                assert(Vertex < 4);
                
                //compile time if to get the shape function I want (just returning barycentric coordinates for this
                //case
                static_if<(Vertex > 0)>([&](auto f){
                    return m_T.row(Vertex-1)*(Eigen::Map<Eigen::Matrix<DataType, 3,1> >(x) - m_x3);
                }).else_([&](auto f) {
                    return 1.0 - (m_T*(Eigen::Map<Eigen::Matrix<DataType, 3,1> >(x) - m_x3)).sum();
                });
                
                
            }
            
            template<unsigned int Vertex>
            inline std::array<DataType, 3> dphi(double *x) {
                
                std::array<DataType, 3> temp;
                static_if<(Vertex >0)>([&](auto f){
                    temp[0] = m_T.row(Vertex-1)(0);
                    temp[1] = m_T.row(Vertex-1)(1);
                    temp[2] = m_T.row(Vertex-1)(2);
                }).else_([&](auto f) {
                    temp[0] = (-m_T.row(0)-m_T.row(1)-m_T.row(2))(0);
                    temp[1] = (-m_T.row(0)-m_T.row(1)-m_T.row(2))(1);
                    temp[2] = (-m_T.row(0)-m_T.row(1)-m_T.row(2))(2);
                });
                
                return temp;
            }
            
        
            //Kinematics stuff (break into new template class then it suits me
            template<typename Matrix>
            inline void F(Matrix &output, State<DataType> &state) {
                
                Eigen::Map<Eigen::VectorXd> m_q0(mapDOFEigen(*m_qDofs[0], state));
                Eigen::Map<Eigen::VectorXd> m_q1(mapDOFEigen(*m_qDofs[1], state));
                Eigen::Map<Eigen::VectorXd> m_q2(mapDOFEigen(*m_qDofs[2], state));
                Eigen::Map<Eigen::VectorXd> m_q3(mapDOFEigen(*m_qDofs[3], state));
                output = m_q0*Eigen::Map<Eigen::VectorXd>(dphi<0>).transpose()+
                         m_q1*Eigen::Map<Eigen::VectorXd>(dphi<1>).transpose()+
                         m_q2*Eigen::Map<Eigen::VectorXd>(dphi<2>).transpose()+
                         m_q3*Eigen::Map<Eigen::VectorXd>(dphi<3>).transpose();
                
            }
        
        protected:
            
            Eigen::Matrix<DataType, 3,3> m_T;
            Eigen::Matrix<DataType, 3,1> m_x3; //maybe don't save this just replace with a pointer and a map
            
            std::array<DOFBase<DataType,0> *, 4> m_qDofs;
            std::array<DOFBase<DataType,1> *, 4> m_qDotDofs;
        
        private:
            
        };
    }
}


#endif /* ShapeFunction_h */
