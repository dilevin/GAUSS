namespace Gauss {
    namespace FEM {
        
        //build a little object that computes neohookean using principal stretches
        struct PSNeohookean {
            
            PSNeohookean() {
                double youngsModulus = 2e6;
                double poissonsRatio = 0.45;
                m_D = 0.5*(youngsModulus*poissonsRatio)/((1.0+poissonsRatio)*(1.0-2.0*poissonsRatio));
                m_C = 0.5*youngsModulus/(2.0*(1.0+poissonsRatio));
            }
            
            template<typename DataType>
            inline void setParameters(DataType youngsModulus, DataType poissonsRatio) {
                m_D = 0.5*(youngsModulus*poissonsRatio)/((1.0+poissonsRatio)*(1.0-2.0*poissonsRatio));
                m_C = 0.5*youngsModulus/(2.0*(1.0+poissonsRatio));
            }
            
            template<typename Derived>
            inline typename Derived::Scalar energy(const Eigen::MatrixBase<Derived> &stretches) {
                
                typename Derived::Scalar J = stretches(0)*stretches(1)*stretches(2);
                
                if(std::abs(J) < 1e-6)
                {
                    J = 1e-3;
                }
                
                typename Derived::Scalar I = stretches(0)*stretches(0) + stretches(1)*stretches(1) + stretches(2)*stretches(2);
                typename Derived::Scalar energy = m_C*((1.0/stablePow(J, static_cast<typename Derived::Scalar>(2.0)))*I - 3.0) + m_D*(J-1.0)*(J-1.0);
                
                return energy;
            }
            
            template<typename Derived>
            inline Eigen::Vector3x<typename Derived::Scalar> gradient(const Eigen::MatrixBase<Derived> &stretches) {
                
                Eigen::Vector3x<typename Derived::Scalar> tmp;
                typename Derived::Scalar I = stretches(0)*stretches(0) + stretches(1)*stretches(1) + stretches(2)*stretches(2);
                typename Derived::Scalar J = stretches(0)*stretches(1)*stretches(2);
                
                if(std::abs(J) < 1e-6)
                {
                    J = 1e-3;
                }
                typename Derived::Scalar J53 = 1.0/stablePow(J,static_cast<typename Derived::Scalar>(5.0));
                typename Derived::Scalar J23 = 1.0/stablePow(J,static_cast<typename Derived::Scalar>(2.0));
                
                tmp[0] = -(2.0/3.0)*m_C*J53*stretches(1)*stretches(2)*I + 2.0*m_C*J23*stretches(0) + 2.0*m_D*(J-1.0)*stretches(1)*stretches(2);
                tmp[1] = -(2.0/3.0)*m_C*J53*stretches(0)*stretches(2)*I + 2.0*m_C*J23*stretches(1) + 2.0*m_D*(J-1.0)*stretches(0)*stretches(2);
                tmp[2] = -(2.0/3.0)*m_C*J53*stretches(0)*stretches(1)*I + 2.0*m_C*J23*stretches(2) + 2.0*m_D*(J-1.0)*stretches(0)*stretches(1);
                
                return tmp;
            }
            
            template<typename Derived>
            inline Eigen::Matrix33x<typename Derived::Scalar> hessian(Eigen::MatrixBase<Derived> &stretches) {
                
                Eigen::Matrix33x<typename Derived::Scalar> tmp;
                
                using DataType = typename Derived::Scalar;
                
                DataType s0 = stretches(0);
                DataType s1 = stretches(1);
                DataType s2 = stretches(2);
                DataType C = m_C;
                DataType D = m_D;
                
                tmp(0,0) = 2*D*std::pow(s1,static_cast<DataType>(2.0))*std::pow(s2,static_cast<DataType>(2.0)) - (0.6666666666666665*C)/stablePow(s0*s1*s2, static_cast<DataType>(2.0)) + (1.111111111111111*C*std::pow(s1,static_cast<DataType>(2.0))*std::pow(s2,static_cast<DataType>(2.0))*(std::pow(s0, static_cast<DataType>(2.0)) + std::pow(s1,static_cast<DataType>(2.0)) + std::pow(s2,static_cast<DataType>(2.0))))/stablePow(s0*s1*s2, static_cast<DataType>(8.0));
                
                tmp(1,0) = (s2*(-0.8888888888888891*C*std::pow(s0, static_cast<DataType>(2.0)) - 0.8888888888888891*C*std::pow(s1,static_cast<DataType>(2.0)) + 0.4444444444444443*C*std::pow(s2,static_cast<DataType>(2.0)) - 2.*D*stablePow(s0*s1*s2, static_cast<DataType>(5.0)) + 4.*D*stablePow(s0*s1*s2, static_cast<DataType>(8.0))))/stablePow(s0*s1*s2, static_cast<DataType>(5.0));
                
                tmp(2,0) = (s1*(-0.8888888888888891*C*std::pow(s0, static_cast<DataType>(2.0)) + 0.4444444444444443*C*std::pow(s1,static_cast<DataType>(2.0)) - 0.8888888888888891*C*std::pow(s2,static_cast<DataType>(2.0)) - 2.*D*stablePow(s0*s1*s2, static_cast<DataType>(5.0)) + 4.*D*stablePow(s0*s1*s2, static_cast<DataType>(8.0))))/stablePow(s0*s1*s2, static_cast<DataType>(5.0));
                
                tmp(0,1) = (s2*(-0.8888888888888891*C*std::pow(s0, static_cast<DataType>(2.0)) - 0.8888888888888891*C*std::pow(s1,static_cast<DataType>(2.0)) + 0.4444444444444443*C*std::pow(s2,static_cast<DataType>(2.0)) - 2.*D*stablePow(s0*s1*s2, static_cast<DataType>(5.0)) + 4.*D*stablePow(s0*s1*s2, static_cast<DataType>(8.0))))/stablePow(s0*s1*s2, static_cast<DataType>(5.0));
                
                tmp(1,1) = 2*D*std::pow(s0, static_cast<DataType>(2.0))*std::pow(s2,static_cast<DataType>(2.0)) - (0.6666666666666665*C)/stablePow(s0*s1*s2, static_cast<DataType>(2.0)) + (1.111111111111111*C*std::pow(s0, static_cast<DataType>(2.0))*std::pow(s2,static_cast<DataType>(2.0))*(std::pow(s0, static_cast<DataType>(2.0)) + std::pow(s1,static_cast<DataType>(2.0)) + std::pow(s2,static_cast<DataType>(2.0))))/stablePow(s0*s1*s2, static_cast<DataType>(8.0));
                
                tmp(2,1) = (s0*(0.4444444444444443*C*std::pow(s0, static_cast<DataType>(2.0)) - 0.8888888888888891*C*std::pow(s1,static_cast<DataType>(2.0)) - 0.8888888888888891*C*std::pow(s2,static_cast<DataType>(2.0)) - 2.*D*stablePow(s0*s1*s2, static_cast<DataType>(5.0)) + 4.*D*stablePow(s0*s1*s2, static_cast<DataType>(8.0))))/stablePow(s0*s1*s2, static_cast<DataType>(5.0));
                
                tmp(0,2) = (s1*(-0.8888888888888891*C*std::pow(s0, static_cast<DataType>(2.0)) + 0.4444444444444443*C*std::pow(s1,static_cast<DataType>(2.0)) - 0.8888888888888891*C*std::pow(s2,static_cast<DataType>(2.0)) - 2.*D*stablePow(s0*s1*s2, static_cast<DataType>(5.0)) + 4.*D*stablePow(s0*s1*s2, static_cast<DataType>(8.0))))/stablePow(s0*s1*s2, static_cast<DataType>(5.0));
                
                tmp(1,2) = (s0*(0.4444444444444443*C*std::pow(s0, static_cast<DataType>(2.0)) - 0.8888888888888891*C*std::pow(s1,static_cast<DataType>(2.0)) - 0.8888888888888891*C*std::pow(s2,static_cast<DataType>(2.0)) - 2.*D*stablePow(s0*s1*s2, static_cast<DataType>(5.0)) + 4.*D*stablePow(s0*s1*s2, static_cast<DataType>(8.0))))/stablePow(s0*s1*s2, static_cast<DataType>(5.0));
                
                tmp(2,2) = 2*D*std::pow(s0, static_cast<DataType>(2.0))*std::pow(s1,static_cast<DataType>(2.0)) - (0.6666666666666665*C)/stablePow(s0*s1*s2, static_cast<DataType>(2.0)) + (1.111111111111111*C*std::pow(s0, static_cast<DataType>(2.0))*std::pow(s1,static_cast<DataType>(2.0))*(std::pow(s0, static_cast<DataType>(2.0)) + std::pow(s1,static_cast<DataType>(2.0)) + std::pow(s2,static_cast<DataType>(2.0))))/stablePow(s0*s1*s2, static_cast<DataType>(8.0));
                
                return tmp;
                
            }
            
            double m_C, m_D;
            
        };
        
        //build a little object that computes neohookean using principal stretches
        struct PSARAP {
            
            PSARAP() {
                m_k = 1e7;
            }
            
            template<typename DataType>
            inline void setParameters(DataType k) {
                m_k = k;
            }
            
            template<typename Derived>
            inline typename Derived::Scalar energy(const Eigen::MatrixBase<Derived> &S) {
                
                typename Derived::Scalar energy = m_k * ((S[0] - static_cast<typename Derived::Scalar>(1)) * (S[0] - static_cast<typename Derived::Scalar>(1)) +
                                                         (S[1] - static_cast<typename Derived::Scalar>(1)) * (S[1] - static_cast<typename Derived::Scalar>(1)) +
                                                         (S[2] - static_cast<typename Derived::Scalar>(1)) * (S[2] - static_cast<typename Derived::Scalar>(1)));
                
                
                return energy;
            }
            
            template<typename Derived>
            inline Eigen::Vector3x<typename Derived::Scalar> gradient(const Eigen::MatrixBase<Derived> &S) {
                
                Eigen::Vector3x<typename Derived::Scalar> tmp;
               
                tmp[0] = static_cast<typename Derived::Scalar>(2)*m_k * (S[0] - static_cast<typename Derived::Scalar>(1));
                tmp[1] = static_cast<typename Derived::Scalar>(2)*m_k * (S[1] - static_cast<typename Derived::Scalar>(1));
                tmp[2] = static_cast<typename Derived::Scalar>(2)*m_k * (S[2] - static_cast<typename Derived::Scalar>(1));
                
                return tmp;
            }
            
            template<typename Derived>
            inline Eigen::Matrix33x<typename Derived::Scalar> hessian(Eigen::MatrixBase<Derived> &S) {
                
                Eigen::Matrix33x<typename Derived::Scalar> tmp;
            
                tmp.setIdentity();
                tmp *= static_cast<typename Derived::Scalar>(2)*m_k;
                
                return tmp;
                
            }
            
            double m_k;
            
        };
        
        /*struct PSCorotatedLinear {
            
            PSCorotatedLinear() {
                double youngsModulus = 2e6;
                double poissonsRatio = 0.45;
                m_lambda = (youngsModulus*poissonsRatio)/((1.0+poissonsRatio)*(1.0-2.0*poissonsRatio));
                m_mu = 0.5*youngsModulus/(2.0*(1.0+poissonsRatio));
                
            }
            
            template<typename DataType>
            inline void setParameters(DataType  youngsModulus, DataType poissonsRatio) {
                m_lambda = (youngsModulus*poissonsRatio)/((1.0+poissonsRatio)*(1.0-2.0*poissonsRatio));
                m_mu = 0.5*youngsModulus/(2.0*(1.0+poissonsRatio));
            }
            
            template<typename Derived>
            inline typename Derived::Scalar energy(const Eigen::MatrixBase<Derived> &S) {
                
                using  Scalar = typename Derived::Scalar;
                
                Scalar energy = m_mu * ((S[0] - 1) * (S[0] - 1) + (S[1] - 1) * (S[1] - 1) + (S[2] - 1) * (S[2] - 1) ) + 0.5*m_lambda *( S[0] + S[1] + S[2] - 3) *( S[0] + S[1] + S[2] - 3.0) ;
                
                return energy;
            }
            
            template<typename Derived>
            inline Eigen::Vector3x<typename Derived::Scalar> gradient(const Eigen::MatrixBase<Derived> &S) {
                
                Eigen::Vector3x<typename Derived::Scalar> tmp;
                
                tmp[0] = static_cast<typename Derived::Scalar>(2) * (S[0] - static_cast<typename Derived::Scalar>(1));
                tmp[1] = static_cast<typename Derived::Scalar>(2) * (S[1] - static_cast<typename Derived::Scalar>(1));
                tmp[2] = static_cast<typename Derived::Scalar>(2) * (S[2] - static_cast<typename Derived::Scalar>(1));
                
                return tmp;
            }
            
            template<typename Derived>
            inline Eigen::Matrix33x<typename Derived::Scalar> hessian(Eigen::MatrixBase<Derived> &S) {
                
                Eigen::Matrix33x<typename Derived::Scalar> tmp;
                
                tmp.setIdentity();
                
                
                return tmp;
                
            }
            
            double m_mu, m_lambda;
        };*/DOFS
    }
}
