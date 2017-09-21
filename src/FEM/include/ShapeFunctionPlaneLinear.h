//
//  ShapeFunctionPlaneLinear.h
//  Gauss
//
//  Created by David Levin on 9/19/17.
//
//

#ifndef ShapeFunctionPlaneLinear_h
#define ShapeFunctionPlaneLinear_h

#include <array>
#include <Utilities.h>

//Plane strain element (implemented as a triangle with a defined thickness)
namespace Gauss {
    namespace FEM {
        
        //Linear Tetrahedral Shape Function
        template<typename DataType>
        class ShapeFunctionPlaneLinear {
        public:
            
            using MatrixB = Eigen::Matrix<DataType, 6,9>;
            using MatrixJ = Eigen::Matrix<DataType, 3, 9>;
            using VectorQ = Eigen::Matrix<DataType, 9,1>;
            
            //convenience matrix type
            template<unsigned int Rows>
            using MatrixDOF = Eigen::Matrix<DataType, Rows, 9>;
            
            template<typename QDOFList, typename QDotDOFList>
            ShapeFunctionPlaneLinear(Eigen::MatrixXd &V, Eigen::MatrixXi &F, QDOFList &qDOFList, QDotDOFList &qDotDOFList) : m_T(2,2) {
                
                if(F.cols() != 3) {
                    std::cout<<"Number of face/vertex indices != 3: Plane strain elements are triangular\n";
                    exit(0);
                }
                
                //currently I assume triangle is in the x,y plane, z direction is just set to be thic
                m_thickness = 1; //set to 1m for now
                
                //build up stuff I need for barycentric coordinates
                m_T <<  (Vert(1,0) - Vert(0,0)), (Vert(2,0) - Vert(0,0)),
                (Vert(1,1) - Vert(0,1)), (Vert(2,1) - Vert(0,1));
                
                
                //invert
                m_x1 << Vert(1,0),
                        Vert(1,1),
                        Vert(1,2);
                m_x2 << Vert(2,0),
                        Vert(2,1),
                        Vert(2,2);
                m_x3 << Vert(0,0),
                        Vert(0,1),
                        Vert(0,2);
                
                m_T = m_T.inverse().eval();
                
                //for the time being assume things come as a stack (qdofs and qdotdofs)
                m_qDofs[0] = qDOFList[0];
                m_qDofs[1] = qDOFList[1];
                m_qDofs[2] = qDOFList[2];
                
                m_qDotDofs[0] = qDotDOFList[0];
                m_qDotDofs[1] = qDotDOFList[1];
                m_qDotDofs[2] = qDotDOFList[2];
            }
            
            //drop my gets for this just because my hands get tired of typing it all the time
            template<unsigned int Vertex>
            inline double phi(double *x) {
                assert(Vertex >= 0);
                assert(Vertex < 3);
                
                DataType phiOut;
                
                //compile time if to get the shape function I want (just returning barycentric coordinates for this
                //case
                static_if<(Vertex > 0)>([&](auto f){
                    phiOut = m_T.row(Vertex-1)*((Eigen::Map<Eigen::Matrix<DataType, 3,1> >(x) - m_x3).head(2));
                }).else_([&](auto f) {
                    phiOut = 1.0 - (m_T*((Eigen::Map<Eigen::Matrix<DataType, 3,1> >(x) - m_x3).head(2))).sum();
                });
                
                return phiOut;
            }
            
            template<unsigned int Vertex>
            inline std::array<DataType, 3> dphi(double *x) {
                
                std::array<DataType, 3> temp;
                static_if<(Vertex >0)>([&](auto f){
                    temp[0] = m_T.row(Vertex-1)(0);
                    temp[1] = m_T.row(Vertex-1)(1);
                    temp[2] = 0.0;
                }).else_([&](auto f) {
                    temp[0] = (-m_T.row(0)-m_T.row(1))(0);
                    temp[1] = (-m_T.row(0)-m_T.row(1))(1);
                    temp[2] = 0;
                });
                
                return temp;
            }
            
            //convert from x in range 0 ->1 to position in element
            inline Eigen::Vector3x<DataType> x(double alphaX, double alphaY, double alphaZ) const {
                return m_x3 + alphaX*(m_x1-m_x3) +alphaX*alphaY*(m_x2-m_x1);
            }
            
            //Kinematics stuff (break into new template class then it suits me
            template<typename Matrix>
            inline void F(Matrix &output, State<DataType> &state) {
                
                Eigen::Map<Eigen::VectorXd> m_q0(mapDOFEigen(*m_qDofs[0], state));
                Eigen::Map<Eigen::VectorXd> m_q1(mapDOFEigen(*m_qDofs[1], state));
                Eigen::Map<Eigen::VectorXd> m_q2(mapDOFEigen(*m_qDofs[2], state));
                output =    m_q0*Eigen::Map<Eigen::VectorXd>(dphi<0>).transpose()+
                            m_q1*Eigen::Map<Eigen::VectorXd>(dphi<1>).transpose()+
                            m_q2*Eigen::Map<Eigen::VectorXd>(dphi<2>).transpose();
                
            }
            
            //Jacobian: derivative with respect to degrees of freedom
            //Jacobian: derivative with respect to degrees of freedom
            inline MatrixJ J(double *x, const State<DataType> &state) {
                
                MatrixJ output;
                
                //just a 3x12 matrix of shape functions
                //kind of assuming everything is initialized before we get here
                double phi0 = phi<0>(x);
                double phi1 = phi<1>(x);
                double phi2 = phi<2>(x);
                
                output.resize(3,9);
                output.setZero();
                output.block(0,0, 3,3) = phi0*Eigen::Matrix<DataType,3,3>::Identity();
                output.block(0,3, 3,3) = phi1*Eigen::Matrix<DataType,3,3>::Identity();
                output.block(0,6, 3,3) = phi2*Eigen::Matrix<DataType,3,3>::Identity();
            
                return output;
                
            }
            
            //spatial derivative of Jacobian (needed to evaluate spatial gradient of displacement/positions
            //rows are derivatives in different directions
            inline MatrixJ GradJ(unsigned int component, double *x, const State<DataType> &state) {
                
                MatrixJ tmp;
                
                tmp.setZero();
                tmp.col(component) = Eigen::Map3x<DataType>(dphi<0>(x).data());
                tmp.col(3+component) = Eigen::Map3x<DataType>(dphi<1>(x).data());
                tmp.col(6+component) = Eigen::Map3x<DataType>(dphi<2>(x).data());
                
                return tmp;
                
            }

            inline VectorQ q(const State<DataType> &state) {
                
                VectorQ tmp;
                
                tmp<<  mapDOFEigen(*m_qDofs[0], state),
                mapDOFEigen(*m_qDofs[1], state),
                mapDOFEigen(*m_qDofs[2], state);
                
                return tmp;
            }

            
            inline double volume() { return 0.5*m_T.determinant()*m_thickness; }
            
            constexpr unsigned int getNumVerts() { return 3; }
            
        protected:
            
            double m_thickness;
            Eigen::Matrix<DataType, 2,2> m_T;
            Eigen::Matrix<DataType, 3,1> m_x1, m_x2, m_x3; //maybe don't save this just replace with a pointer and a map
            
            std::array<DOFBase<DataType,0> *, 3> m_qDofs;
            std::array<DOFBase<DataType,1> *, 3> m_qDotDofs;
            
        private:
            
        };
        
    }
}

#endif /* ShapeFunctionPlaneLinear_h */
