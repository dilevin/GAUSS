//
//  Energy.h
//  Gauss
//
//  Created by David Levin on 3/13/17.
//
//

#ifndef Energy_h
#define Energy_h

namespace Gauss {
    namespace FEM {
        
        /////// Kinetic Energy Terms //////
        template<typename DataType, typename ShapeFunction>
        class EnergyKineticNone : public ShapeFunction {
    
        };
        
        //Non-lumped kinetic energy
        template<typename DataType, typename ShapeFunction>
        class EnergyKineticNonLumped : public ShapeFunction {
        public:
            template<typename QDOFList, typename QDotDOFList>
            EnergyKineticNonLumped(Eigen::MatrixXd &V, Eigen::MatrixXi &F, QDOFList &qDOFList, QDotDOFList &qDotDOFList) : ShapeFunction(V, F, qDOFList, qDotDOFList) {
                m_rho = 1000.0;
            }
            
            inline void setDensity(double density) { m_rho = density; }
            
            inline double getDensity() { return m_rho; }
            
            inline DataType getValue(double *x, const State<DataType> &state) {
                auto v = ShapeFunction::J(x, state)*ShapeFunction::qDot();
                return 0.5*v.transpose()*v;
            }
            
            //infinitessimal gradients and hessians
            template<typename Vector>
            inline void getGradient(Vector &f, double *x, const State<DataType> &state) {
                //do Nothing (I don't think I use this for anything really)
            }
            
            template<typename Matrix>
            inline void getHessian(Matrix &H, double *x, const State<DataType> &state) {
             
                //compute mass matrix at point X = m_rho*J(x)^T*J(x)
                typename ShapeFunction::MatrixJ M = ShapeFunction::J(x, state);
                    
                H = m_rho*M.transpose()*M;
            }
            
        protected:
            double m_rho;
            
        private:
            
        };
        
        /////// Potential Energy Terms //////
        template<typename DataType, typename ShapeFunction>
        class EnergyPotentialNone : public ShapeFunction {
        public:
            template<typename DOFList>
            EnergyPotentialNone(Eigen::MatrixXd &V, Eigen::MatrixXi &F, DOFList &dofList) : ShapeFunction(V,F, dofList) { }
        };
        
        template<typename DataType, typename ShapeFunction>
        class EnergyLinearElasticity : public ShapeFunction {
        public:
            template<typename QDOFList, typename QDotDOFList>
            EnergyLinearElasticity(Eigen::MatrixXd &V, Eigen::MatrixXi &F, QDOFList &qDOFList, QDotDOFList &qDotDOFList) : ShapeFunction(V, F, qDOFList, qDotDOFList) {
                setParameters(1e6, 0.45);
                
            }
            
            inline void setParameters(double youngsModulus, double poissonsRatio) {
                m_E = youngsModulus;
                m_mu = poissonsRatio;
                
                m_C.setZero();
                m_C(0,0) = 1.0-m_mu;
                m_C(0,1) = m_mu;
                m_C(0,2) = m_mu;
                m_C(1,0) = m_mu;
                m_C(1,1) = 1.0-m_mu;
                m_C(1,2) = m_mu;
                m_C(2,0) = m_mu;
                m_C(2,1) = m_mu;
                m_C(2,2) = 1.0-m_mu;
                m_C(3,3) = 0.5*(1.0-2.0*m_mu);
                m_C(4,4) = 0.5*(1.0-2.0*m_mu);
                m_C(5,5) = 0.5*(1.0-2.0*m_mu);
                m_C *= (m_E/((1.0+m_mu)*(1.0-2.0*m_mu)));
            }
            
            inline DataType getValue(double *x, const State<DataType> &state) {
                
                auto q = ShapeFunction::q();
                
                return -0.5*q.transpose()*B(this, x, state).transpose()*m_C*B(this, x, state)*q;
            }
            
            template<typename Vector>
            inline void getGradient(Vector &f, double *x, const State<DataType> &state) {
                
                /*//returning the force which is really the negative gradient*/
                f = -B(this, x, state).transpose()*m_C*B(this, x, state)*ShapeFunction::q(state);

            }
            
            template<typename Matrix>
            inline void getHessian(Matrix &H, double *x, const State<DataType> &state) {
                //H = -ShapeFunction::B(x,state).transpose()*m_C*ShapeFunction::B(x,state);
                
                H = -B(this, x, state).transpose()*m_C*B(this, x, state);
            
            }
            
            template<typename Matrix>
            inline void getCauchyStress(Matrix &S, double *x, State<DataType> &state) {
                S = m_C*B(this,x,state)*ShapeFunction::q(state);
            }
            
            inline const DataType & getE() const { return m_E; }
            inline const DataType & getMu() const { return m_mu; }
            
        protected:
            DataType m_E, m_mu;
            Eigen::Matrix66x<DataType> m_C;
            
        private:
            
        };
        
        //////// Body Force Terms /////////
        template<typename DataType, typename ShapeFunction>
        class BodyForceNone : public ShapeFunction {
        public:
             template<typename QDOFList, typename QDotDOFList, typename ...Params>
            BodyForceNone(Eigen::MatrixXd &V, Eigen::MatrixXi &F, QDOFList &qDOFList, QDotDOFList &qDotDOFList, Params ...params) : ShapeFunction(V,F, qDOFList, qDotDOFList) { }
        };
        
        template<typename DataType, typename ShapeFunction>
        class BodyForceGravity : public ShapeFunction {
        public:
            template<typename QDOFList, typename QDotDOFList>
            BodyForceGravity(Eigen::MatrixXd &V, Eigen::MatrixXi &F, QDOFList &qDOFList, QDotDOFList &qDotDOFList) : ShapeFunction(V,F, qDOFList, qDotDOFList)
            {
                m_rho = 1;
                m_g << 0.0,0.0,0.0;
            }
            
            void setParams(double rho, double gx, double gy, double gz) {
                m_rho = rho;
                m_g << gx,gy,gz;
            }
            
            inline  double getValue(double *x, const State<DataType> &state) {
            
                return  -ShapeFunction::q(state).transpose()*ShapeFunction::J(x,state).transpose()*m_rho*m_g;
            }
            
            template<typename Vector>
            inline void getGradient(Vector &f, double *x, const State<DataType> &state) {
            
                //assuming force does same rate of work so generalized force = J'T*f
                f = ShapeFunction::J(x,state).transpose()*m_rho*m_g;
                
            }
            
            template<typename Matrix>
            inline void getHessian(Matrix &H, double *x, const State<DataType> &state) {
                
            }
            
        protected:
            Eigen::Vector3x<DataType> m_g;
            double m_rho;
        };
        
    }
}


#endif /* Energy_h */
