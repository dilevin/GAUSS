#ifndef ENERGY_LOUBIGNAC
#define ENERGY_LOUBIGNAC

#include<cmath>

template<typename DataType, typename ShapeFunction>
class EnergyLoubignac : public virtual ShapeFunction {
public:
    //Only works for linear tetrahedral elements for now
    template<typename QDOFList, typename QDotDOFList>
    EnergyLoubignac(Eigen::MatrixXd &V, Eigen::MatrixXi &F, QDOFList &qDOFList, QDotDOFList &qDotDOFList) : ShapeFunction(V, F, qDOFList, qDotDOFList) {
        
    }

    inline DataType getValue(double *x, const State<DataType> &state) {
        

        std::cout<<"Error: EnergyLoubignac does not actually compute an energy \n";
        return 0.0;
    }
    
    template<typename Vector>
    inline void getGradient(Vector &f, double *x, const State<DataType> &state) {
        
        //RHS for Loubignac (hack to get Gauss to integrate this correctly
        //Acts on mesh using 6x1 stress DOF
        std::tuple<DataType *, unsigned int> ptr = state.template getStatePtr<0>();
        Eigen::VectorXx<DataType> s0 = Eigen::Map<Eigen::VectorXx<DataType> >(&std::get<0>(ptr)[6*(this->q()[0]->getLocalId())/(this->q()[0]->getNumScalarDOF())], 6);
        Eigen::VectorXx<DataType> s1 = Eigen::Map<Eigen::VectorXx<DataType> >(&std::get<0>(ptr)[6*(this->q()[1]->getLocalId())/(this->q()[1]->getNumScalarDOF())], 6);
        Eigen::VectorXx<DataType> s2 = Eigen::Map<Eigen::VectorXx<DataType> >(&std::get<0>(ptr)[6*(this->q()[2]->getLocalId())/(this->q()[2]->getNumScalarDOF())], 6);
        Eigen::VectorXx<DataType> s3 = Eigen::Map<Eigen::VectorXx<DataType> >(&std::get<0>(ptr)[6*(this->q()[3]->getLocalId())/(this->q()[3]->getNumScalarDOF())], 6);
        
        
        Eigen::VectorXx<DataType> fInt;
        fInt.resize(12,1);
        
        fInt =-B(this, x, state).transpose()*(this->template phi<0>(x)*s0+
                                            this->template phi<1>(x)*s1 +
                                            this->template phi<2>(x)*s2 +
                                            this->template phi<3>(x)*s3);
        f = fInt;
    }
    
    template<typename Matrix>
    inline void getHessian(Matrix &H, double *x, const State<DataType> &state) {
        std::cout<<"Error: EnergyLoubignac does not compute a gradient \n";
    }
    
    template<typename Matrix>
    inline void getCauchyStress(Matrix &S, double *x, State<DataType> &state) {
        
    }
    
protected:

private:
    
};

#endif

