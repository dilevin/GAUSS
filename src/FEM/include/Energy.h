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
            
            //infinitessimal gradients and hessians
            template<typename Vector>
            inline void getGradient(Vector &f, double *x, const State<DataType> &state) {
                //do Nothing (I don't think I use this for anything really)
            }
            
            template<typename Matrix>
            inline void getHessian(Matrix &H, double *x, const State<DataType> &state) {
             
                //compute mass matrix at point X = m_rho*J(x)^T*J(x)
                Eigen::MatrixXd M;
                ShapeFunction::J(M, x, state);

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
        
        //Linear Elasticity
        template<typename DataType, typename ShapeFunction>
        class EnergyLinearElasticity : public ShapeFunction {
        public:
            template<typename QDOFList, typename QDotDOFList>
            EnergyLinearElasticity(Eigen::MatrixXd &V, Eigen::MatrixXi &F, QDOFList &qDOFList, QDotDOFList &qDotDOFList) : ShapeFunction(V, F, qDOFList, qDotDOFList) {
                m_E = 1e7;
                m_mu = 0.45;
            }
            
            inline void setParameters(double youngsModulus, double poissonsRatio) {
                m_E = youngsModulus;
                m_mu = poissonsRatio;
            }
            
            template<typename Vector>
            inline void getGradient(Vector &f, double *x, const State<DataType> &state) {
                
            }
            
            template<typename Matrix>
            inline void getHessian(Matrix &H, double *x, const State<DataType> &state) {
                
            }
            
            //accessors
            inline DataType & setE() { return m_E; }
            inline DataType & setMu() { return m_mu; }
            
            inline const DataType & getE() const { return m_E; }
            inline const DataType & getMu() const { return m_mu; }
            
        protected:
            DataType m_E, m_mu;
            
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
            
            template<typename Vector>
            inline void getGradient(Vector &f, double *x, const State<DataType> &state) {
                f = m_rho*m_g;
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
