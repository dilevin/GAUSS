#ifndef ENERGY_PRINCIPALSTRETCH
#define ENERGY_PRINCIPALSTRETCH

#include<cmath>
#include <igl/svd3x3.h>
//Energy PS is an object that can compute the energy, gradient and hessian of the energy, expressed as a function of principal stretches.
template<typename DataType, typename ShapeFunction, typename EnergyPS>
class EnergyPrincipalStretch : public virtual ShapeFunction {
public:
    template<typename QDOFList, typename QDotDOFList>
    EnergyPrincipalStretch(Eigen::MatrixXd &V, Eigen::MatrixXi &F, QDOFList &qDOFList, QDotDOFList &qDotDOFList) : ShapeFunction(V, F, qDOFList, qDotDOFList) {
        
        
    }

    inline EnergyPS & getPrincipalStretchObject() { return m_ps; }
    
    inline DataType getValue(double *x, const State<DataType> &state) {
        
        Eigen::Matrix33x<float> F = (ShapeFunction::F(x,state) + Eigen::Matrix<DataType,3,3>::Identity()).template cast<float>();;
        
        igl::svd3x3(F,m_U,m_S,m_V);
        
        //evaluate energy
        return static_cast<DataType>(m_ps.energy(m_S));
    }
    
    template<typename Vector>
    inline void getGradient(Vector &f, double *x, const State<DataType> &state) {
        
        //Deformation Gradient
        Eigen::Matrix33x<float> F = (ShapeFunction::F(x,state) + Eigen::Matrix<DataType,3,3>::Identity()).template cast<float>();
        
        //Eigen::JacobiSVD<Eigen::MatrixXx<DataType> > svd;
        //DataType maxVal = F.maxCoeff();
        //svd.compute(F/maxVal, Eigen::ComputeFullU | Eigen::ComputeFullV);
        
        igl::svd3x3(F,m_U,m_S,m_V);
        
        //Eigen::Vector3x<DataType> Plam = m_ps.gradient(svd.singularValues()*maxVal);
        Eigen::Vector3x<float> Plam = m_ps.gradient(m_S);
        
        //Eigen::MatrixXx<DataType> P = svd.matrixU()*(Plam.asDiagonal()*svd.matrixV().transpose());
        Eigen::MatrixXx<float> P = m_U*(Plam.asDiagonal()*m_V.transpose());
        P.transposeInPlace();
        P.resize(9,1);
                        
        
        //build force vector
        f = -ShapeFunction::GradJ(0,x,state).transpose()*P.block(0,0, 3,1).template cast<DataType>() +
            -ShapeFunction::GradJ(1,x,state).transpose()*P.block(3,0, 3,1).template cast<DataType>() +
            -ShapeFunction::GradJ(2,x,state).transpose()*P.block(6,0, 3,1).template cast<DataType>();
        
        //Deformation Gradient
        /*double f11, f12, f13, f21, f22, f23, f31, f32, f33;
        
        Eigen::Matrix<DataType, 3,3> F2 = ShapeFunction::F(x,state);
        
        f11 = F2(0,0)+1.0;
        f12 = F2(0,1);
        f13 = F2(0,2);
        f21 = F2(1,0);
        f22 = F2(1,1)+1.0;
        f23 = F2(1,2);
        f31 = F2(2,0);
        f32 = F2(2,1);
        f33 = F2(2,2)+1.0;
        
        //Force Vector computation from Mathematica notebook
        Eigen::VectorXx<DataType> f2;
        Eigen::Matrix<DataType, 9,1> dw;
        double C,D;
        
        double youngsModulus = 2e6;
        double poissonsRatio = 0.45;
        D = 0.5*(youngsModulus*poissonsRatio)/((1.0+poissonsRatio)*(1.0-2.0*poissonsRatio));
        C = 0.5*youngsModulus/(2.0*(1.0+poissonsRatio));
        
        dw[0] = 2*D*(f23*f32 - f22*f33)*(1 + f13*f22*f31 - f12*f23*f31 - f13*f21*f32 + f11*f23*f32 + f12*f21*f33 - f11*f22*f33) + (2*C*(3*f11*(-(f13*f22*f31) + f12*f23*f31 + f13*f21*f32 - f12*f21*f33) + std::pow(f11,2)*(-2*f23*f32 + 2*f22*f33) + (f23*f32 - f22*f33)*(std::pow(f12,2) + std::pow(f13,2) + std::pow(f21,2) + std::pow(f22,2) + std::pow(f23,2) + std::pow(f31,2) + std::pow(f32,2) + std::pow(f33,2))))/(3.*stablePow(-(f13*f22*f31) + f12*f23*f31 + f13*f21*f32 - f11*f23*f32 - f12*f21*f33 + f11*f22*f33,5.0));
        
        dw[1] = 2*D*(f23*f31 - f21*f33)*(-1 - f13*f22*f31 + f12*f23*f31 + f13*f21*f32 - f11*f23*f32 - f12*f21*f33 + f11*f22*f33) - (2*C*(std::pow(f12,2)*(-2*f23*f31 + 2*f21*f33) + 3*f12*(f13*f22*f31 - f13*f21*f32 + f11*f23*f32 - f11*f22*f33) + (f23*f31 - f21*f33)*(std::pow(f11,2) + std::pow(f13,2) + std::pow(f21,2) + std::pow(f22,2) + std::pow(f23,2) + std::pow(f31,2) + std::pow(f32,2) + std::pow(f33,2))))/(3.*stablePow(-(f13*f22*f31) + f12*f23*f31 + f13*f21*f32 - f11*f23*f32 - f12*f21*f33 + f11*f22*f33,5.0));
        
        dw[2] = 2*D*(f22*f31 - f21*f32)*(1 + f13*f22*f31 - f12*f23*f31 - f13*f21*f32 + f11*f23*f32 + f12*f21*f33 - f11*f22*f33) + (2*C*f13)/stablePow(-(f13*f22*f31) + f12*f23*f31 + f13*f21*f32 - f11*f23*f32 - f12*f21*f33 + f11*f22*f33,2.0) + (2*C*(f22*f31 - f21*f32)*(std::pow(f11,2) + std::pow(f12,2) + std::pow(f13,2) + std::pow(f21,2) + std::pow(f22,2) + std::pow(f23,2) + std::pow(f31,2) + std::pow(f32,2) + std::pow(f33,2)))/(3.*stablePow(-(f13*f22*f31) + f12*f23*f31 + f13*f21*f32 - f11*f23*f32 - f12*f21*f33 + f11*f22*f33,5.0));
        
        dw[3] = 2*D*(f13*f32 - f12*f33)*(-1 - f13*f22*f31 + f12*f23*f31 + f13*f21*f32 - f11*f23*f32 - f12*f21*f33 + f11*f22*f33) + (2*C*f21)/stablePow(-(f13*f22*f31) + f12*f23*f31 + f13*f21*f32 - f11*f23*f32 - f12*f21*f33 + f11*f22*f33,2.0) - (2*C*(f13*f32 - f12*f33)*(std::pow(f11,2) + std::pow(f12,2) + std::pow(f13,2) + std::pow(f21,2) + std::pow(f22,2) + std::pow(f23,2) + std::pow(f31,2) + std::pow(f32,2) + std::pow(f33,2)))/(3.*stablePow(-(f13*f22*f31) + f12*f23*f31 + f13*f21*f32 - f11*f23*f32 - f12*f21*f33 + f11*f22*f33,5.0));
        
        dw[4] = 2*D*(f13*f31 - f11*f33)*(1 + f13*f22*f31 - f12*f23*f31 - f13*f21*f32 + f11*f23*f32 + f12*f21*f33 - f11*f22*f33) + (2*C*f22)/stablePow(-(f13*f22*f31) + f12*f23*f31 + f13*f21*f32 - f11*f23*f32 - f12*f21*f33 + f11*f22*f33,2.0) + (2*C*(f13*f31 - f11*f33)*(std::pow(f11,2) + std::pow(f12,2) + std::pow(f13,2) + std::pow(f21,2) + std::pow(f22,2) + std::pow(f23,2) + std::pow(f31,2) + std::pow(f32,2) + std::pow(f33,2)))/(3.*stablePow(-(f13*f22*f31) + f12*f23*f31 + f13*f21*f32 - f11*f23*f32 - f12*f21*f33 + f11*f22*f33,5.0));
        
        dw[5] = 2*D*(f12*f31 - f11*f32)*(-1 - f13*f22*f31 + f12*f23*f31 + f13*f21*f32 - f11*f23*f32 - f12*f21*f33 + f11*f22*f33) + (2*C*f23)/stablePow(-(f13*f22*f31) + f12*f23*f31 + f13*f21*f32 - f11*f23*f32 - f12*f21*f33 + f11*f22*f33,2.0) - (2*C*(f12*f31 - f11*f32)*(std::pow(f11,2) + std::pow(f12,2) + std::pow(f13,2) + std::pow(f21,2) + std::pow(f22,2) + std::pow(f23,2) + std::pow(f31,2) + std::pow(f32,2) + std::pow(f33,2)))/(3.*stablePow(-(f13*f22*f31) + f12*f23*f31 + f13*f21*f32 - f11*f23*f32 - f12*f21*f33 + f11*f22*f33,5.0));
        
        dw[6] = 2*D*(f13*f22 - f12*f23)*(1 + f13*f22*f31 - f12*f23*f31 - f13*f21*f32 + f11*f23*f32 + f12*f21*f33 - f11*f22*f33) + (2*C*f31)/stablePow(-(f13*f22*f31) + f12*f23*f31 + f13*f21*f32 - f11*f23*f32 - f12*f21*f33 + f11*f22*f33,2.0) + (2*C*(f13*f22 - f12*f23)*(std::pow(f11,2) + std::pow(f12,2) + std::pow(f13,2) + std::pow(f21,2.0) + std::pow(f22,2) + std::pow(f23,2) + std::pow(f31,2) + std::pow(f32,2) + std::pow(f33,2)))/(3.*stablePow(-(f13*f22*f31) + f12*f23*f31 + f13*f21*f32 - f11*f23*f32 - f12*f21*f33 + f11*f22*f33,5.0));
        
        dw[7] = 2*D*(f13*f21 - f11*f23)*(-1 - f13*f22*f31 + f12*f23*f31 + f13*f21*f32 - f11*f23*f32 - f12*f21*f33 + f11*f22*f33) + (2*C*f32)/stablePow(-(f13*f22*f31) + f12*f23*f31 + f13*f21*f32 - f11*f23*f32 - f12*f21*f33 + f11*f22*f33,2.0) - (2*C*(f13*f21 - f11*f23)*(std::pow(f11,2) + std::pow(f12,2) + std::pow(f13,2) + std::pow(f21,2) + std::pow(f22,2) + std::pow(f23,2) + std::pow(f31,2) + std::pow(f32,2) + std::pow(f33,2)))/(3.*stablePow(-(f13*f22*f31) + f12*f23*f31 + f13*f21*f32 - f11*f23*f32 - f12*f21*f33 + f11*f22*f33,5.0));
        
        dw[8] = 2*D*(f12*f21 - f11*f22)*(1 + f13*f22*f31 - f12*f23*f31 - f13*f21*f32 + f11*f23*f32 + f12*f21*f33 - f11*f22*f33) - (2*C*((-(f12*f21) + f11*f22)*(std::pow(f11,2) + std::pow(f12,2) + std::pow(f13,2) + std::pow(f21,2) + std::pow(f22,2) + std::pow(f23,2) + std::pow(f31,2) + std::pow(f32,2)) + 3*(f13*f22*f31 - f12*f23*f31 - f13*f21*f32 + f11*f23*f32)*f33 + 2*(f12*f21 - f11*f22)*std::pow(f33,2)))/(3.*stablePow(-(f13*f22*f31) + f12*f23*f31 + f13*f21*f32 - f11*f23*f32 - f12*f21*f33 + f11*f22*f33,5.0));
        
        f2 = -ShapeFunction::GradJ(0,x,state).transpose()*dw.segment(0,3) +
        -ShapeFunction::GradJ(1,x,state).transpose()*dw.segment(3,3) +
        -ShapeFunction::GradJ(2,x,state).transpose()*dw.segment(6,3);
        
        //P(3) = dw[3];
        //build force vector
        f = -ShapeFunction::GradJ(0,x,state).transpose()*P.block(0,0, 3,1) +
        -ShapeFunction::GradJ(1,x,state).transpose()*P.block(3,0, 3,1) +
        -ShapeFunction::GradJ(2,x,state).transpose()*P.block(6,0, 3,1);
        
        /*std::cout<<"DW: "<<dw.transpose()<<"\n";
        std::cout<<"P: "<<P.transpose()<<"\n";
        std::cout<<"ERROR: "<<(f-f2).norm()<<"\n";*/
        
        
    }
    
    template<typename Matrix>
    inline void getHessian(Matrix &H, double *x, const State<DataType> &state) {
        //H = -ShapeFunction::B(x,state).transpose()*m_C*ShapeFunction::B(x,state);
        /*Eigen::Matrix<DataType, 3,3> F = ShapeFunction::F(x,state);
        double f11, f12, f13, f21, f22, f23, f31, f32, f33;
        
        f11 = F(0,0)+1.0;
        f12 = F(0,1);
        f13 = F(0,2);
        f21 = F(1,0);
        f22 = F(1,1)+1.0;
        f23 = F(1,2);
        f31 = F(2,0);
        f32 = F(2,1);
        f33 = F(2,2)+1.0;
        
        //H = -B(this, x, state).transpose()*m_C*B(this, x, state);
        Eigen::Matrix<DataType,9,9> ddw;
        typename ShapeFunction::MatrixJ gradX, gradY, gradZ;
        gradX = ShapeFunction::GradJ(0,x,state);
        gradY = ShapeFunction::GradJ(1,x,state);
        gradZ = ShapeFunction::GradJ(2,x,state);*/
        
        
        //No hessian yet
        typename ShapeFunction::MatrixJ J;
        H.resize(J.cols(), J.cols());
        H.setZero();
        //End mathematica code
        /*H = -gradX.transpose()*ddw.block(0,0,3,3)*gradX +
        -gradX.transpose()*ddw.block(0,3,3,3)*gradY +
        -gradX.transpose()*ddw.block(0,6,3,3)*gradZ +
        -gradY.transpose()*ddw.block(3,0,3,3)*gradX +
        -gradY.transpose()*ddw.block(3,3,3,3)*gradY +
        -gradY.transpose()*ddw.block(3,6,3,3)*gradZ +
        -gradZ.transpose()*ddw.block(6,0,3,3)*gradX +
        -gradZ.transpose()*ddw.block(6,3,3,3)*gradY +
        -gradZ.transpose()*ddw.block(6,6,3,3)*gradZ;*/
        
        
    }
    
    template<typename Matrix>
    inline void getCauchyStress(Matrix &S, double *x, State<DataType> &state) {
        std::cout<<"GetCauchyStress not implemented in EnergyPrincipalStrain \n";
        std::exit(1);
    }
    
    inline const DataType getE() const { return 10.0; }
    
protected:
    
    //Energy PS object
    EnergyPS m_ps;
    Eigen::Vector3x<float> m_S;
    Eigen::Matrix33x<float> m_U,m_V;
    
private:
    
};

#endif


