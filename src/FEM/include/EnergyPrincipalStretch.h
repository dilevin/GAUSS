#ifndef ENERGY_PRINCIPALSTRETCH
#define ENERGY_PRINCIPALSTRETCH

#include<cmath>
#include <igl/svd3x3.h>
//Energy PS is an object that can compute the energy, gradient and hessian of the energy, expressed as a function of principal stretches.
template<typename DataType, typename ShapeFunction, typename EnergyPS>
class EnergyPrincipalStretch : public virtual ShapeFunction {
public:
    template<typename QDOFList, typename QDotDOFList>
    EnergyPrincipalStretch(Eigen::MatrixXx<DataType> &V, Eigen::MatrixXi &F, QDOFList &qDOFList, QDotDOFList &qDotDOFList) : ShapeFunction(V, F, qDOFList, qDotDOFList) {
        
        
    }

    inline EnergyPS & getPrincipalStretchObject() { return m_ps; }
    
    template<typename ...Params>
    inline void setParameters(Params ...params) {
        m_ps.setParameters(params...);
    }
    
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
        
        igl::svd3x3(F,m_U,m_S,m_V);
        
        //check multiplicity, small random permutation if eigenvalues not unique
        if(std::fabs(m_S[0] - m_S[1]) < 1e-5 || std::fabs(m_S[1] - m_S[2]) < 1e-5 || std::fabs(m_S[0] - m_S[2]) < 1e-5) {
            F += Eigen::Matrix33x<float>::Random()*1e-5;
            igl::svd3x3(F,m_U,m_S,m_V);
        }
        
        Eigen::Vector3x<float> Plam = m_ps.gradient(m_S);
        
        //Eigen::MatrixXx<DataType> P = svd.matrixU()*(Plam.asDiagonal()*svd.matrixV().transpose());
        Eigen::MatrixXx<float> P = m_U*(Plam.asDiagonal()*m_V.transpose());
        P.transposeInPlace();
        P.resize(9,1);
                        
        
        //build force vector
        f = -ShapeFunction::GradJ(0,x,state).transpose()*P.block(0,0, 3,1).template cast<DataType>() +
            -ShapeFunction::GradJ(1,x,state).transpose()*P.block(3,0, 3,1).template cast<DataType>() +
            -ShapeFunction::GradJ(2,x,state).transpose()*P.block(6,0, 3,1).template cast<DataType>();
        
    }
    
    template<typename Matrix>
    inline void getHessian(Matrix &H, DataType *x, const State<DataType> &state) {
        
        //Deformation Gradient
        Eigen::Matrix33x<float> F = (ShapeFunction::F(x,state) + Eigen::Matrix<DataType,3,3>::Identity()).template cast<float>();
        
        igl::svd3x3(F,m_U,m_S,m_V);
        
        //Eigen::Vector3x<DataType> Plam = m_ps.gradient(svd.singularValues()*maxVal);
        Eigen::Matrix<DataType,9,9, Eigen::RowMajor> ddw2;
        Eigen::Vector3x<float> Plam = m_ps.gradient(m_S);
        Eigen::Matrix33x<float> Plam2 = m_ps.hessian(m_S);
        
        //derivative of SVD wrt to F
        Eigen::dSVD(m_dU, m_dS, m_dV, m_U, m_S, m_V);
        Eigen::Matrix33x<DataType> rowMat;
        
        //This formatting is quite complicated
        for(unsigned int r = 0; r <3; ++r) {
            for(unsigned int s = 0; s<3; ++s) {
                Eigen::Vector3x<float> PlamVec = Plam2*m_dS[r][s];
                rowMat = (m_dU[r][s]*Plam.asDiagonal()*m_V.transpose() + m_U*Plam.asDiagonal()*m_dV[r][s].transpose() + m_U*PlamVec.asDiagonal()*m_V.transpose()).template cast<DataType>();
                rowMat.transposeInPlace();
                ddw2.row(3*r + s) = Eigen::Map<Eigen::Matrix<DataType, 1,9> >(rowMat.data(), 9);
                
            }
        }

        typename ShapeFunction::MatrixJ gradX, gradY, gradZ;
        gradX = ShapeFunction::GradJ(0,x,state);
        gradY = ShapeFunction::GradJ(1,x,state);
        gradZ = ShapeFunction::GradJ(2,x,state);
        
        H = -gradX.transpose()*ddw2.block(0,0,3,3)*gradX +
        -gradX.transpose()*ddw2.block(0,3,3,3)*gradY +
        -gradX.transpose()*ddw2.block(0,6,3,3)*gradZ +
        -gradY.transpose()*ddw2.block(3,0,3,3)*gradX +
        -gradY.transpose()*ddw2.block(3,3,3,3)*gradY +
        -gradY.transpose()*ddw2.block(3,6,3,3)*gradZ +
        -gradZ.transpose()*ddw2.block(6,0,3,3)*gradX +
        -gradZ.transpose()*ddw2.block(6,3,3,3)*gradY +
        -gradZ.transpose()*ddw2.block(6,6,3,3)*gradZ;
        
    }
    
    template<typename Matrix>
    inline void getCauchyStress(Matrix &S, double *x, State<DataType> &state) {
        std::cout<<"GetCauchyStress not implemented in EnergyPrincipalStrain \n";
        std::exit(1);
    }
    
    inline const EnergyPS & material() { m_ps; }
    
    inline const DataType getE() const { return 10.0; }
    
protected:
    
    //Energy PS object
    EnergyPS m_ps;
    Eigen::Vector3x<float> m_S;
    Eigen::Matrix33x<float> m_U,m_V;
    Eigen::Tensor3333x<float> m_dU, m_dV;
    Eigen::Tensor333x<float> m_dS;
private:
    
};

#endif


