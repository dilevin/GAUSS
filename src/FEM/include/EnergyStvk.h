#ifndef ENERGY_STVK
#define ENERGY_STVK

#include <cmath>

template<typename DataType, typename ShapeFunction>
class EnergyStvk : public virtual ShapeFunction{
public:
    template<typename QDOFList, typename QDotDOFList>
    EnergyStvk(Eigen::MatrixXd &V, Eigen::MatrixXi &F, QDOFList &qDOFList, QDotDOFList &qDotDOFList) : ShapeFunction(V, F, qDOFList, qDotDOFList) {
        setParameters(2e6, 0.45);
    }

    inline void setParameters(double youngsModulus, double poissonsRatio) {
        m_D = 0.5*(youngsModulus*poissonsRatio)/((1.0+poissonsRatio)*(1.0-2.0*poissonsRatio));
        m_C = 0.5*youngsModulus/(2.0*(1.0+poissonsRatio));
    }

    //Stvk Energy
    inline DataType getValue(double *x, const State<DataType> &state) {
    	// std::cout<<"Stvk Energy"<<std::endl;
        Eigen::Matrix<DataType, 3,3> F = ShapeFunction::F(x,state) + Eigen::Matrix<DataType,3,3>::Identity();
        //E = greenstrain
        Eigen::Matrix<DataType, 3,3> E = 0.5*(F.transpose()*F - Eigen::Matrix<DataType,3,3>::Identity());
        double Psi_energyDensity = 2*m_C*(E*E).trace() + m_D*E.trace()*E.trace();
        return Psi_energyDensity;
    }
    
    template<typename Vector>
    inline void getGradient(Vector &f, double *x, const State<DataType> &state) {
    	// std::cout<<"stvk gradient"<<std::endl;
		Eigen::Matrix<DataType, 3,3> F = ShapeFunction::F(x,state) + Eigen::Matrix<DataType,3,3>::Identity();
		//Deformation Gradient
        double f11, f12, f13, f21, f22, f23, f31, f32, f33;
        f11 = F(0,0);
        f12 = F(0,1);
        f13 = F(0,2);
        f21 = F(1,0);
        f22 = F(1,1);
        f23 = F(1,2);
        f31 = F(2,0);
        f32 = F(2,1);
        f33 = F(2,2);
        // Eigen::Matrix<DataType, 3,3> E = 0.5*(F.transpose()*F) - Eigen::Matrix<DataType,3,3>::Identity();
        // Eigen::Matrix<DataType, 3,3> P = F*(2*m_mu*E + m_lamda*E.trace()*Eigen::Matrix<DataType,3,3>::Identity());
    	// Eigen::Matrix<DataType, 3,3> dEdx = -1*ShapeFunction::volume()*P*ShapeFunction::getInvRefShapeMatrix();
		
		Eigen::Matrix<DataType, 9,1> dw;
		double C,D;
        C = m_C;
        D = m_D;

		dw[0] = 2.*C*f11*(-1. + std::pow(f11, 2) + std::pow(f21, 2) + std::pow(f31, 2)) + 1.*D*f11*(-1. + std::pow(f11, 2) + std::pow(f21, 2) + std::pow(f31, 2)) + 2.*C*f12*(f11*f12 + f21*f22 + f31*f32) + 1.*D*f11*(-1. + std::pow(f12, 2) + std::pow(f22, 2) + std::pow(f32, 2)) + 2.*C*f13*(f11*f13 + f21*f23 + f31*f33) + 1.*D*f11*(-1. + std::pow(f13, 2) + std::pow(f23, 2) + std::pow(f33, 2));
    	dw[1] = 1.*D*f12*(-1. + std::pow(f11, 2) + std::pow(f21, 2) + std::pow(f31, 2)) + 2.*C*f11*(f11*f12 + f21*f22 + f31*f32) + 2.*C*f12*(-1. + std::pow(f12, 2) + std::pow(f22, 2) + std::pow(f32, 2)) + 1.*D*f12*(-1. + std::pow(f12, 2) + std::pow(f22, 2) + std::pow(f32, 2)) + 2.*C*f13*(f12*f13 + f22*f23 + f32*f33) + 1.*D*f12*(-1. + std::pow(f13, 2) + std::pow(f23, 2) + std::pow(f33, 2));
    	dw[2] = 1.*D*f13*(-1. + std::pow(f11, 2) + std::pow(f21, 2) + std::pow(f31, 2)) + 1.*D*f13*(-1. + std::pow(f12, 2) + std::pow(f22, 2) + std::pow(f32, 2)) + 2.*C*f11*(f11*f13 + f21*f23 + f31*f33) + 2.*C*f12*(f12*f13 + f22*f23 + f32*f33) + 2.*C*f13*(-1. + std::pow(f13, 2) + std::pow(f23, 2) + std::pow(f33, 2)) + 1.*D*f13*(-1. + std::pow(f13, 2) + std::pow(f23, 2) + std::pow(f33, 2));
    	dw[3] = 2.*C*f21*(-1. + std::pow(f11, 2) + std::pow(f21, 2) + std::pow(f31, 2)) + 1.*D*f21*(-1. + std::pow(f11, 2) + std::pow(f21, 2) + std::pow(f31, 2)) + 2.*C*f22*(f11*f12 + f21*f22 + f31*f32) + 1.*D*f21*(-1. + std::pow(f12, 2) + std::pow(f22, 2) + std::pow(f32, 2)) + 2.*C*f23*(f11*f13 + f21*f23 + f31*f33) + 1.*D*f21*(-1. + std::pow(f13, 2) + std::pow(f23, 2) + std::pow(f33, 2));
    	dw[4] = 1.*D*f22*(-1. + std::pow(f11, 2) + std::pow(f21, 2) + std::pow(f31, 2)) + 2.*C*f21*(f11*f12 + f21*f22 + f31*f32) + 2.*C*f22*(-1. + std::pow(f12, 2) + std::pow(f22, 2) + std::pow(f32, 2)) + 1.*D*f22*(-1. + std::pow(f12, 2) + std::pow(f22, 2) + std::pow(f32, 2)) + 2.*C*f23*(f12*f13 + f22*f23 + f32*f33) + 1.*D*f22*(-1. + std::pow(f13, 2) + std::pow(f23, 2) + std::pow(f33, 2));
    	dw[5] = 1.*D*f23*(-1. + std::pow(f11, 2) + std::pow(f21, 2) + std::pow(f31, 2)) + 1.*D*f23*(-1. + std::pow(f12, 2) + std::pow(f22, 2) + std::pow(f32, 2)) + 2.*C*f21*(f11*f13 + f21*f23 + f31*f33) + 2.*C*f22*(f12*f13 + f22*f23 + f32*f33) + 2.*C*f23*(-1. + std::pow(f13, 2) + std::pow(f23, 2) + std::pow(f33, 2)) + 1.*D*f23*(-1. + std::pow(f13, 2) + std::pow(f23, 2) + std::pow(f33, 2));
    	dw[6] = 2.*C*f31*(-1. + std::pow(f11, 2) + std::pow(f21, 2) + std::pow(f31, 2)) + 1.*D*f31*(-1. + std::pow(f11, 2) + std::pow(f21, 2) + std::pow(f31, 2)) + 2.*C*f32*(f11*f12 + f21*f22 + f31*f32) + 1.*D*f31*(-1. + std::pow(f12, 2) + std::pow(f22, 2) + std::pow(f32, 2)) + 2.*C*f33*(f11*f13 + f21*f23 + f31*f33) + 1.*D*f31*(-1. + std::pow(f13, 2) + std::pow(f23, 2) + std::pow(f33, 2));
    	dw[7] = 1.*D*(-1. + std::pow(f11, 2) + std::pow(f21, 2) + std::pow(f31, 2))*f32 + 2.*C*f31*(f11*f12 + f21*f22 + f31*f32) + 2.*C*f32*(-1. + std::pow(f12, 2) + std::pow(f22, 2) + std::pow(f32, 2)) + 1.*D*f32*(-1. + std::pow(f12, 2) + std::pow(f22, 2) + std::pow(f32, 2)) + 2.*C*f33*(f12*f13 + f22*f23 + f32*f33) + 1.*D*f32*(-1. + std::pow(f13, 2) + std::pow(f23, 2) + std::pow(f33, 2));
    	dw[8] = 1.*D*(-1. + std::pow(f11, 2) + std::pow(f21, 2) + std::pow(f31, 2))*f33 + 1.*D*(-1. + std::pow(f12, 2) + std::pow(f22, 2) + std::pow(f32, 2))*f33 + 2.*C*f31*(f11*f13 + f21*f23 + f31*f33) + 2.*C*f32*(f12*f13 + f22*f23 + f32*f33) + 2.*C*f33*(-1. + std::pow(f13, 2) + std::pow(f23, 2) + std::pow(f33, 2)) + 1.*D*f33*(-1. + std::pow(f13, 2) + std::pow(f23, 2) + std::pow(f33, 2));

    	f = -ShapeFunction::GradJ(0,x,state).transpose()*dw.segment(0,3) +
            -ShapeFunction::GradJ(1,x,state).transpose()*dw.segment(3,3) +
            -ShapeFunction::GradJ(2,x,state).transpose()*dw.segment(6,3);
    }
    
    template<typename Matrix>
    inline void getHessian(Matrix &H, double *x, const State<DataType> &state) {
    	// std::cout<<"Stvk Hessian"<<std::endl;
    	Eigen::Matrix<DataType, 3,3> F = ShapeFunction::F(x,state) + Eigen::Matrix<DataType,3,3>::Identity();
		//Deformation Gradient
        double f11, f12, f13, f21, f22, f23, f31, f32, f33;
        f11 = F(0,0);
        f12 = F(0,1);
        f13 = F(0,2);
        f21 = F(1,0);
        f22 = F(1,1);
        f23 = F(1,2);
        f31 = F(2,0);
        f32 = F(2,1);
        f33 = F(2,2);

        Eigen::Matrix<DataType,9,9> ddw;
        typename ShapeFunction::MatrixJ gradX, gradY, gradZ;
        gradX = ShapeFunction::GradJ(0,x,state);
        gradY = ShapeFunction::GradJ(1,x,state);
        gradZ = ShapeFunction::GradJ(2,x,state);
        
        double C,D;
        C = m_C;
        D = m_D;

    	ddw(0,0) = 4.*C*std::pow(f11, 2) + 2.*D*std::pow(f11, 2) + 2.*C*std::pow(f12, 2) + 2.*C*std::pow(f13, 2) + 2.*C*(-1. + std::pow(f11, 2) + std::pow(f21, 2) + std::pow(f31, 2)) + 1.*D*(-1. + std::pow(f11, 2) + std::pow(f21, 2) + std::pow(f31, 2)) + 1.*D*(-1. + std::pow(f12, 2) + std::pow(f22, 2) + std::pow(f32, 2)) + 1.*D*(-1. + std::pow(f13, 2) + std::pow(f23, 2) + std::pow(f33, 2));
		ddw(1,0) = 2.*D*f11*f12 + C*(4.*f11*f12 + 2.*f21*f22 + 2.*f31*f32);
		ddw(2,0) = 2.*D*f11*f13 + C*(4.*f11*f13 + 2.*f21*f23 + 2.*f31*f33);
		ddw(3,0) = 2.*D*f11*f21 + C*(4.*f11*f21 + 2.*f12*f22 + 2.*f13*f23);
		ddw(4,0) = 2.*C*f12*f21 + 2.*D*f11*f22;
		ddw(5,0) = 2.*C*f13*f21 + 2.*D*f11*f23;
		ddw(6,0) = 2.*D*f11*f31 + C*(4.*f11*f31 + 2.*f12*f32 + 2.*f13*f33);
		ddw(7,0) = 2.*C*f12*f31 + 2.*D*f11*f32;
		ddw(8,0) = 2.*C*f13*f31 + 2.*D*f11*f33;
		ddw(0,1) = 2.*D*f11*f12 + C*(4.*f11*f12 + 2.*f21*f22 + 2.*f31*f32);
		ddw(1,1) = 2.*C*std::pow(f11, 2) + 4.*C*std::pow(f12, 2) + 2.*D*std::pow(f12, 2) + 2.*C*std::pow(f13, 2) + 1.*D*(-1. + std::pow(f11, 2) + std::pow(f21, 2) + std::pow(f31, 2)) + 2.*C*(-1. + std::pow(f12, 2) + std::pow(f22, 2) + std::pow(f32, 2)) + 1.*D*(-1. + std::pow(f12, 2) + std::pow(f22, 2) + std::pow(f32, 2)) + 1.*D*(-1. + std::pow(f13, 2) + std::pow(f23, 2) + std::pow(f33, 2));
		ddw(2,1) = 2.*D*f12*f13 + C*(4.*f12*f13 + 2.*f22*f23 + 2.*f32*f33);
		ddw(3,1) = 2.*D*f12*f21 + 2.*C*f11*f22;
		ddw(4,1) = 2.*D*f12*f22 + C*(2.*f11*f21 + 4.*f12*f22 + 2.*f13*f23);
		ddw(5,1) = 2.*C*f13*f22 + 2.*D*f12*f23;
		ddw(6,1) = 2.*D*f12*f31 + 2.*C*f11*f32;
		ddw(7,1) = 2.*D*f12*f32 + C*(2.*f11*f31 + 4.*f12*f32 + 2.*f13*f33);
		ddw(8,1) = 2.*C*f13*f32 + 2.*D*f12*f33;
		ddw(0,2) = 2.*D*f11*f13 + C*(4.*f11*f13 + 2.*f21*f23 + 2.*f31*f33);
		ddw(1,2) = 2.*D*f12*f13 + C*(4.*f12*f13 + 2.*f22*f23 + 2.*f32*f33);
		ddw(2,2) = 2.*C*std::pow(f11, 2) + 2.*C*std::pow(f12, 2) + 4.*C*std::pow(f13, 2) + 2.*D*std::pow(f13, 2) + 1.*D*(-1. + std::pow(f11, 2) + std::pow(f21, 2) + std::pow(f31, 2)) + 1.*D*(-1. + std::pow(f12, 2) + std::pow(f22, 2) + std::pow(f32, 2)) + 2.*C*(-1. + std::pow(f13, 2) + std::pow(f23, 2) + std::pow(f33, 2)) + 1.*D*(-1. + std::pow(f13, 2) + std::pow(f23, 2) + std::pow(f33, 2));
		ddw(3,2) = 2.*D*f13*f21 + 2.*C*f11*f23;
		ddw(4,2) = 2.*D*f13*f22 + 2.*C*f12*f23;
		ddw(5,2) = 2.*D*f13*f23 + C*(2.*f11*f21 + 2.*f12*f22 + 4.*f13*f23);
		ddw(6,2) = 2.*D*f13*f31 + 2.*C*f11*f33;
		ddw(7,2) = 2.*D*f13*f32 + 2.*C*f12*f33;
		ddw(8,2) = 2.*D*f13*f33 + C*(2.*f11*f31 + 2.*f12*f32 + 4.*f13*f33);
		ddw(0,3) = 2.*D*f11*f21 + C*(4.*f11*f21 + 2.*f12*f22 + 2.*f13*f23);
		ddw(1,3) = 2.*D*f12*f21 + 2.*C*f11*f22;
		ddw(2,3) = 2.*D*f13*f21 + 2.*C*f11*f23;
		ddw(3,3) = 4.*C*std::pow(f21, 2) + 2.*D*std::pow(f21, 2) + 2.*C*std::pow(f22, 2) + 2.*C*std::pow(f23, 2) + 2.*C*(-1. + std::pow(f11, 2) + std::pow(f21, 2) + std::pow(f31, 2)) + 1.*D*(-1. + std::pow(f11, 2) + std::pow(f21, 2) + std::pow(f31, 2)) + 1.*D*(-1. + std::pow(f12, 2) + std::pow(f22, 2) + std::pow(f32, 2)) + 1.*D*(-1. + std::pow(f13, 2) + std::pow(f23, 2) + std::pow(f33, 2));
		ddw(4,3) = 2.*D*f21*f22 + C*(2.*f11*f12 + 4.*f21*f22 + 2.*f31*f32);
		ddw(5,3) = 2.*D*f21*f23 + C*(2.*f11*f13 + 4.*f21*f23 + 2.*f31*f33);
		ddw(6,3) = 2.*D*f21*f31 + C*(4.*f21*f31 + 2.*f22*f32 + 2.*f23*f33);
		ddw(7,3) = 2.*C*f22*f31 + 2.*D*f21*f32;
		ddw(8,3) = 2.*C*f23*f31 + 2.*D*f21*f33;
		ddw(0,4) = 2.*C*f12*f21 + 2.*D*f11*f22;
		ddw(1,4) = 2.*D*f12*f22 + C*(2.*f11*f21 + 4.*f12*f22 + 2.*f13*f23);
		ddw(2,4) = 2.*D*f13*f22 + 2.*C*f12*f23;
		ddw(3,4) = 2.*D*f21*f22 + C*(2.*f11*f12 + 4.*f21*f22 + 2.*f31*f32);
		ddw(4,4) = 2.*C*std::pow(f21, 2) + 4.*C*std::pow(f22, 2) + 2.*D*std::pow(f22, 2) + 2.*C*std::pow(f23, 2) + 1.*D*(-1. + std::pow(f11, 2) + std::pow(f21, 2) + std::pow(f31, 2)) + 2.*C*(-1. + std::pow(f12, 2) + std::pow(f22, 2) + std::pow(f32, 2)) + 1.*D*(-1. + std::pow(f12, 2) + std::pow(f22, 2) + std::pow(f32, 2)) + 1.*D*(-1. + std::pow(f13, 2) + std::pow(f23, 2) + std::pow(f33, 2));
		ddw(5,4) = 2.*D*f22*f23 + C*(2.*f12*f13 + 4.*f22*f23 + 2.*f32*f33);
		ddw(6,4) = 2.*D*f22*f31 + 2.*C*f21*f32;
		ddw(7,4) = 2.*D*f22*f32 + C*(2.*f21*f31 + 4.*f22*f32 + 2.*f23*f33);
		ddw(8,4) = 2.*C*f23*f32 + 2.*D*f22*f33;
		ddw(0,5) = 2.*C*f13*f21 + 2.*D*f11*f23;
		ddw(1,5) = 2.*C*f13*f22 + 2.*D*f12*f23;
		ddw(2,5) = 2.*D*f13*f23 + C*(2.*f11*f21 + 2.*f12*f22 + 4.*f13*f23);
		ddw(3,5) = 2.*D*f21*f23 + C*(2.*f11*f13 + 4.*f21*f23 + 2.*f31*f33);
		ddw(4,5) = 2.*D*f22*f23 + C*(2.*f12*f13 + 4.*f22*f23 + 2.*f32*f33);
		ddw(5,5) = 2.*C*std::pow(f21, 2) + 2.*C*std::pow(f22, 2) + 4.*C*std::pow(f23, 2) + 2.*D*std::pow(f23, 2) + 1.*D*(-1. + std::pow(f11, 2) + std::pow(f21, 2) + std::pow(f31, 2)) + 1.*D*(-1. + std::pow(f12, 2) + std::pow(f22, 2) + std::pow(f32, 2)) + 2.*C*(-1. + std::pow(f13, 2) + std::pow(f23, 2) + std::pow(f33, 2)) + 1.*D*(-1. + std::pow(f13, 2) + std::pow(f23, 2) + std::pow(f33, 2));
		ddw(6,5) = 2.*D*f23*f31 + 2.*C*f21*f33;
		ddw(7,5) = 2.*D*f23*f32 + 2.*C*f22*f33;
		ddw(8,5) = 2.*D*f23*f33 + C*(2.*f21*f31 + 2.*f22*f32 + 4.*f23*f33);
		ddw(0,6) = 2.*D*f11*f31 + C*(4.*f11*f31 + 2.*f12*f32 + 2.*f13*f33);
		ddw(1,6) = 2.*D*f12*f31 + 2.*C*f11*f32;
		ddw(2,6) = 2.*D*f13*f31 + 2.*C*f11*f33;
		ddw(3,6) = 2.*D*f21*f31 + C*(4.*f21*f31 + 2.*f22*f32 + 2.*f23*f33);
		ddw(4,6) = 2.*D*f22*f31 + 2.*C*f21*f32;
		ddw(5,6) = 2.*D*f23*f31 + 2.*C*f21*f33;
		ddw(6,6) = 4.*C*std::pow(f31, 2) + 2.*D*std::pow(f31, 2) + 2.*C*(-1. + std::pow(f11, 2) + std::pow(f21, 2) + std::pow(f31, 2)) + 1.*D*(-1. + std::pow(f11, 2) + std::pow(f21, 2) + std::pow(f31, 2)) + 2.*C*std::pow(f32, 2) + 1.*D*(-1. + std::pow(f12, 2) + std::pow(f22, 2) + std::pow(f32, 2)) + 2.*C*std::pow(f33, 2) + 1.*D*(-1. + std::pow(f13, 2) + std::pow(f23, 2) + std::pow(f33, 2));
		ddw(7,6) = 2.*D*f31*f32 + C*(2.*f11*f12 + 2.*f21*f22 + 4.*f31*f32);
		ddw(8,6) = 2.*D*f31*f33 + C*(2.*f11*f13 + 2.*f21*f23 + 4.*f31*f33);
		ddw(0,7) = 2.*C*f12*f31 + 2.*D*f11*f32;
		ddw(1,7) = 2.*D*f12*f32 + C*(2.*f11*f31 + 4.*f12*f32 + 2.*f13*f33);
		ddw(2,7) = 2.*D*f13*f32 + 2.*C*f12*f33;
		ddw(3,7) = 2.*C*f22*f31 + 2.*D*f21*f32;
		ddw(4,7) = 2.*D*f22*f32 + C*(2.*f21*f31 + 4.*f22*f32 + 2.*f23*f33);
		ddw(5,7) = 2.*D*f23*f32 + 2.*C*f22*f33;
		ddw(6,7) = 2.*D*f31*f32 + C*(2.*f11*f12 + 2.*f21*f22 + 4.*f31*f32);
		ddw(7,7) = 2.*C*std::pow(f31, 2) + 1.*D*(-1. + std::pow(f11, 2) + std::pow(f21, 2) + std::pow(f31, 2)) + 4.*C*std::pow(f32, 2) + 2.*D*std::pow(f32, 2) + 2.*C*(-1. + std::pow(f12, 2) + std::pow(f22, 2) + std::pow(f32, 2)) + 1.*D*(-1. + std::pow(f12, 2) + std::pow(f22, 2) + std::pow(f32, 2)) + 2.*C*std::pow(f33, 2) + 1.*D*(-1. + std::pow(f13, 2) + std::pow(f23, 2) + std::pow(f33, 2));
		ddw(8,7) = 2.*D*f32*f33 + C*(2.*f12*f13 + 2.*f22*f23 + 4.*f32*f33);
		ddw(0,8) = 2.*C*f13*f31 + 2.*D*f11*f33;
		ddw(1,8) = 2.*C*f13*f32 + 2.*D*f12*f33;
		ddw(2,8) = 2.*D*f13*f33 + C*(2.*f11*f31 + 2.*f12*f32 + 4.*f13*f33);
		ddw(3,8) = 2.*C*f23*f31 + 2.*D*f21*f33;
		ddw(4,8) = 2.*C*f23*f32 + 2.*D*f22*f33;
		ddw(5,8) = 2.*D*f23*f33 + C*(2.*f21*f31 + 2.*f22*f32 + 4.*f23*f33);
		ddw(6,8) = 2.*D*f31*f33 + C*(2.*f11*f13 + 2.*f21*f23 + 4.*f31*f33);
		ddw(7,8) = 2.*D*f32*f33 + C*(2.*f12*f13 + 2.*f22*f23 + 4.*f32*f33);
		ddw(8,8) = 2.*C*std::pow(f31, 2) + 1.*D*(-1. + std::pow(f11, 2) + std::pow(f21, 2) + std::pow(f31, 2)) + 2.*C*std::pow(f32, 2) + 1.*D*(-1. + std::pow(f12, 2) + std::pow(f22, 2) + std::pow(f32, 2)) + 4.*C*std::pow(f33, 2) + 2.*D*std::pow(f33, 2) + 2.*C*(-1. + std::pow(f13, 2) + std::pow(f23, 2) + std::pow(f33, 2)) + 1.*D*(-1. + std::pow(f13, 2) + std::pow(f23, 2) + std::pow(f33, 2));
    
		//End mathematica code
        H = -gradX.transpose()*ddw.block(0,0,3,3)*gradX +
            -gradX.transpose()*ddw.block(0,3,3,3)*gradY +
            -gradX.transpose()*ddw.block(0,6,3,3)*gradZ +
            -gradY.transpose()*ddw.block(3,0,3,3)*gradX +
            -gradY.transpose()*ddw.block(3,3,3,3)*gradY +
            -gradY.transpose()*ddw.block(3,6,3,3)*gradZ +
            -gradZ.transpose()*ddw.block(6,0,3,3)*gradX +
            -gradZ.transpose()*ddw.block(6,3,3,3)*gradY +
            -gradZ.transpose()*ddw.block(6,6,3,3)*gradZ;
    }

	template<typename Matrix>
    inline void getCauchyStress(Matrix &S, double *x, State<DataType> &state) {
        double mu = 2*m_C;
        double lambda = 2*m_D;

        Eigen::Matrix<DataType, 3,3> F = ShapeFunction::F(x,state) + Eigen::Matrix<DataType,3,3>::Identity();
        double J = F.determinant();
        Eigen::Matrix<DataType, 3,3> E = 0.5*(F.transpose()*F) - Eigen::Matrix<DataType,3,3>::Identity();
        Eigen::Matrix<DataType, 3,3> P = F*(m_C*E + 2*m_D*E.trace()*Eigen::Matrix<DataType,3,3>::Identity());

        S = (P * F.transpose()) / J;
    }  


    inline const DataType & getC() const { return m_C; }
    inline const DataType & getD() const { return m_D; }
    
    inline const DataType & getE() const { return m_C; }

protected:
    DataType m_C, m_D;
private:
    
};

#endif
