//
//  ShapeFunctionHexTrilinear.h
//  Gauss
//
//  Created by David Levin on 4/27/17.
//
//

#ifndef ShapeFunctionHexTrilinear_h
#define ShapeFunctionHexTrilinear_h
#include <UtilitiesEigen.h>



namespace Gauss {
    namespace FEM {

        
        //TODO Add trilinear hexahedral shape function + approriate quadrature rules
        template<typename DataType>
        class ShapeFunctionHexTrilinear {
        public:
            
            using MatrixB = Eigen::Matrix<DataType, 6,24>;
            using MatrixJ = Eigen::Matrix<DataType, 3, 24>;
            using VectorQ = Eigen::Matrix<DataType, 24,1>;
            
            //convenience matrix type
            template<unsigned int Rows>
            using MatrixDOF = Eigen::Matrix<DataType, Rows, 24>;
            
            template<typename QDOFList, typename QDotDOFList>
            ShapeFunctionHexTrilinear(Eigen::MatrixXd &V, Eigen::MatrixXi &F, QDOFList &qDOFList, QDotDOFList &qDotDOFList) {

                //for the time being assume things come as a stack (qdofs and qdotdofs)
                m_qDofs[0] = qDOFList[0];
                m_qDofs[1] = qDOFList[1];
                m_qDofs[2] = qDOFList[2];
                m_qDofs[3] = qDOFList[3];
                m_qDofs[4] = qDOFList[4];
                m_qDofs[5] = qDOFList[5];
                m_qDofs[6] = qDOFList[6];
                m_qDofs[7] = qDOFList[7];
                
                
                m_qDotDofs[0] = qDotDOFList[0];
                m_qDotDofs[1] = qDotDOFList[1];
                m_qDotDofs[2] = qDotDOFList[2];
                m_qDotDofs[3] = qDotDOFList[3];
                m_qDotDofs[4] = qDotDOFList[4];
                m_qDotDofs[5] = qDotDOFList[5];
                m_qDotDofs[6] = qDotDOFList[6];
                m_qDotDofs[7] = qDotDOFList[7];
                
                m_x0 << Vert(0,0), Vert(0,1), Vert(0,2);
                m_dx = (V.block(F(6), 0, 1,3) - V.block(F(0), 0, 1, 3)).transpose();
                
                //using node numberings from http://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.Ch11.d/AFEM.Ch11.pdf
                //but with right handed coordinate system
                //     4----------5
                //    /|         /|
                //   7-|--------6 |
                //   | |        | |
                //   | 0 -------|-1    +y
                //   |/         |/     |
                //   3----------2      0-- +x
                //                    /
                //                   +z
                
            }
            
            //drop my gets for this just because my hands get tired of typing it all the time
            template<int Vertex>
            inline double phi(double *x) {
                
                assert(Vertex < 8);
                Eigen::Vector3d e = alpha(x);
                double phiOut = 0.0;
                
                static_if<(Vertex==0)>([&](auto f) {
                    phiOut =  (1.0/8.0)*(1-e(0))*(1-e(1))*(1-e(2));
                });
                static_if<(Vertex==1)>([&](auto f) {
                    phiOut =  (1.0/8.0)*(1+e(0))*(1-e(1))*(1-e(2));
                });
                static_if<(Vertex==2)>([&](auto f) {
                    phiOut =  (1.0/8.0)*(1+e(0))*(1-e(1))*(1+e(2));
                });
                static_if<(Vertex==3)>([&](auto f) {
                    phiOut =  (1.0/8.0)*(1-e(0))*(1-e(1))*(1+e(2));
                });
                static_if<(Vertex==4)>([&](auto f) {
                    phiOut =  (1.0/8.0)*(1-e(0))*(1+e(1))*(1-e(2));
                });
                static_if<Vertex==5>([&](auto f) {
                    phiOut =  (1.0/8.0)*(1+e(0))*(1+e(1))*(1-e(2));
                });
                static_if<(Vertex==6)>([&](auto f) {
                    phiOut =  (1.0/8.0)*(1+e(0))*(1+e(1))*(1+e(2));
                });
                static_if<(Vertex==7)>([&](auto f) {
                    phiOut =  (1.0/8.0)*(1-e(0))*(1+e(1))*(1+e(2));
                });
                
                
                /*static_if<(Vertex==0)>([&](auto f) {
                    return (1.0/8.0)*(1-f(e)(0))*(1-f(e)(1))*(1-f(e)(2));
                }).else_([&](auto f) {
                    static_if<(Vertex==1)>([&](auto f) {
                        return (1.0/8.0)*(1+f(e)(0))*(1-f(e)(1))*(1-f(e)(2));
                    }).else_([&](auto f) {
                        static_if<(Vertex==2)>([&](auto f) {
                            return (1.0/8.0)*(1+e(0))*(1+e(1))*(1-e(2));
                        }).else_([&](auto f) {
                            static_if<(Vertex==3)>([&](auto f) {
                                return (1.0/8.0)*(1-e(0))*(1+e(1))*(1-e(2));
                            }).else_([&](auto f) {
                                static_if<(Vertex==4)>([&](auto f) {
                                    return (1.0/8.0)*(1-e(0))*(1-e(1))*(1+e(2));
                                }).else_([&](auto f) {
                                    static_if<(Vertex==5)>([&](auto f) {
                                        return (1.0/8.0)*(1+f(e)(0))*(1-f(e)(1))*(1+f(e)(2));
                                    }).else_([&](auto f) {
                                        static_if<(Vertex==6)>([&](auto f) {
                                            return (1.0/8.0)*(1+e(0))*(1+e(1))*(1+e(2));
                                        }).else_([&](auto f) {
                                            //final shape function
                                            return (1.0/8.0)*(1-e(0))*(1+e(1))*(1+e(2));
                                        });
                                    });
                                });
                            });
                        });
                    });
                    
                });*/
            
                return phiOut;
            }
            
            inline Eigen::Vector3x<DataType> x(double alphaX, double alphaY, double alphaZ) const {
                return m_x0 + Eigen::Vector3x<DataType>((alphaX+1.)*m_dx(0)*0.5, (alphaY+1.)*m_dx(1)*0.5,(alphaZ+1.)*m_dx(2)*0.5);
            }
            
            inline Eigen::Vector3d alpha(double *x) {
            
                Eigen::Vector3d e = Eigen::Map<Eigen::Vector3d>(x)-m_x0;
                e(0) = 2.0*(e(0)/m_dx(0))-1.0;
                e(1) = 2.0*(e(1)/m_dx(1))-1.0;
                e(2) = 2.0*(e(2)/m_dx(2))-1.0;
                
                return e;
            }
            
            template<unsigned int Vertex>
            inline std::array<DataType, 3> dphi(double *x) {
                
                std::array<DataType,3> deriv;
                
                assert(Vertex < 8);
                Eigen::Vector3d e = alpha(x);
              
                
                //need dalpha/dx
                
                static_if<(Vertex==0)>([&](auto f) {
                    deriv[0] = -0.25*((1. - e(1))*(1. - e(2)))/m_dx(0);
                    deriv[1] = -0.25*((1. - e(0))*(1. - e(2)))/m_dx(1);
                    deriv[2] = -0.25*((1. - e(0))*(1. - e(1)))/m_dx(2);
                });
                static_if<(Vertex==1)>([&](auto f) {
                    deriv[0] =  0.25*((1. - e(1))*(1. - e(2)))/m_dx(0);
                    deriv[1] = -0.25*((1. + e(0))*(1. - e(2)))/m_dx(1);
                    deriv[2] = -0.25*((1. + e(0))*(1. - e(1)))/m_dx(2);
                });
                static_if<(Vertex==2)>([&](auto f) {
                    deriv[0] =  0.25*((1. - e(1))*(1. + e(2)))/m_dx(0);
                    deriv[1] = -0.25*((1. + e(0))*(1. + e(2)))/m_dx(1);
                    deriv[2] =  0.25*((1. + e(0))*(1. - e(1)))/m_dx(2);
                });
                static_if<(Vertex==3)>([&](auto f) {
                    deriv[0] = -0.25*((1. - e(1))*(1. + e(2)))/m_dx(0);
                    deriv[1] = -0.25*((1. - e(0))*(1. + e(2)))/m_dx(1);
                    deriv[2] =  0.25*((1. - e(0))*(1. - e(1)))/m_dx(2);
                });
                static_if<(Vertex==4)>([&](auto f) {
                    deriv[0] = -0.25*((1. + e(1))*(1. - e(2)))/m_dx(0);
                    deriv[1] =  0.25*((1. - e(0))*(1. - e(2)))/m_dx(1);
                    deriv[2] = -0.25*((1. - e(0))*(1. + e(1)))/m_dx(2);
                });
                static_if<Vertex==5>([&](auto f) {
                    deriv[0] =  0.25*((1. + e(1))*(1. - e(2)))/m_dx(0);
                    deriv[1] =  0.25*((1. + e(0))*(1. - e(2)))/m_dx(1);
                    deriv[2] = -0.25*((1. + e(0))*(1. + e(1)))/m_dx(2);
                });
                static_if<(Vertex==6)>([&](auto f) {
                    deriv[0] =  0.25*((1. + e(1))*(1. + e(2)))/m_dx(0);
                    deriv[1] =  0.25*((1. + e(0))*(1. + e(2)))/m_dx(1);
                    deriv[2] =  0.25*((1. + e(0))*(1. + e(1)))/m_dx(2);
                });
                static_if<(Vertex==7)>([&](auto f) {
                    deriv[0] = -0.25*((1. + e(1))*(1. + e(2)))/m_dx(0);
                    deriv[1] =  0.25*((1. - e(0))*(1. + e(2)))/m_dx(1);
                    deriv[2] =  0.25*((1. - e(0))*(1. + e(1)))/m_dx(2);
                });

                
                return deriv;
            }
            
            
            //Kinematics stuff (break into new template class then it suits me
            template<typename Matrix>
            inline void F(Matrix &output, State<DataType> &state) {
                
                std::cout<<"Not implemented yet \n";
                assert(1 == 0);
            }
            
    
            /*inline Eigen::Matrix<DataType, 6, 24> B(double *x, const State<DataType> &state) {
            
                MatrixB tmp;
                
                tmp <<      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
                            0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
                            0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
                            0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
                            0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
                            0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0;
                
                return tmp;
            }*/
            
            inline VectorQ q(const State<DataType> &state) {
                
                VectorQ tmp;
                
                
                tmp<<  mapDOFEigen(*m_qDofs[0], state),
                       mapDOFEigen(*m_qDofs[1], state),
                       mapDOFEigen(*m_qDofs[2], state),
                       mapDOFEigen(*m_qDofs[3], state),
                       mapDOFEigen(*m_qDofs[4], state),
                       mapDOFEigen(*m_qDofs[5], state),
                       mapDOFEigen(*m_qDofs[6], state),
                       mapDOFEigen(*m_qDofs[7], state);
                
                return tmp;
            }
            
            //Jacobian: derivative with respect to degrees of freedom
            inline MatrixJ J(double *x, const State<DataType> &state) {
                
                MatrixJ output;
                
                //just a 3x12 matrix of shape functions
                //kind of assuming everything is initialized before we get here
                double phi0 = phi<0>(x);
                double phi1 = phi<1>(x);
                double phi2 = phi<2>(x);
                double phi3 = phi<3>(x);
                double phi4 = phi<4>(x);
                double phi5 = phi<5>(x);
                double phi6 = phi<6>(x);
                double phi7 = phi<7>(x);

                output.resize(3,24);
                output.setZero();
                output.block(0,0, 3,3) = phi0*Eigen::Matrix<DataType,3,3>::Identity();
                output.block(0,3, 3,3) = phi1*Eigen::Matrix<DataType,3,3>::Identity();
                output.block(0,6, 3,3) = phi2*Eigen::Matrix<DataType,3,3>::Identity();
                output.block(0,9, 3,3) = phi3*Eigen::Matrix<DataType,3,3>::Identity();
                output.block(0,12, 3,3) = phi4*Eigen::Matrix<DataType,3,3>::Identity();
                output.block(0,15, 3,3) = phi5*Eigen::Matrix<DataType,3,3>::Identity();
                output.block(0,18, 3,3) = phi6*Eigen::Matrix<DataType,3,3>::Identity();
                output.block(0,21, 3,3) = phi7*Eigen::Matrix<DataType,3,3>::Identity();
                
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
                tmp.col(9+component) = Eigen::Map3x<DataType>(dphi<3>(x).data());
                tmp.col(12+component) = Eigen::Map3x<DataType>(dphi<4>(x).data());
                tmp.col(15+component) = Eigen::Map3x<DataType>(dphi<5>(x).data());
                tmp.col(18+component) = Eigen::Map3x<DataType>(dphi<6>(x).data());
                tmp.col(21+component) = Eigen::Map3x<DataType>(dphi<7>(x).data());
                
                return tmp;
                
            }
            
            inline double volume() { return m_dx(0)*m_dx(1)*m_dx(2); }
            
            constexpr unsigned int getNumVerts() { return 8; }
            
        protected:
            
            std::array<DOFBase<DataType,0> *, 8> m_qDofs;
            std::array<DOFBase<DataType,1> *, 8> m_qDotDofs;
            Eigen::Vector3d m_dx;
            Eigen::Vector3d m_x0;
            
        private:
            
        };
    }
}

#endif /* ShapeFunctionHexTrilinear_h */
