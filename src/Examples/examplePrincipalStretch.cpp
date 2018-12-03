//Example using principal stretch-based materials
#include <functional>

#include <Qt3DIncludes.h>
#include <GaussIncludes.h>
#include <ParticleSystemIncludes.h>
#include <FEMIncludes.h>
#include <Newton.h>

#include <iostream>
#include <tuple>
//Any extra things I need such as constraints
#include <ConstraintFixedPoint.h>
#include <TimeStepperEulerImplicitBFGS.h>
#include <TimeStepperEulerImplicitBFGS.h>
#include <type_traits>

using namespace Gauss;
using namespace FEM;
using namespace ParticleSystem; //For Force Spring

//build a little object that computes neohookean using principal stretches
struct NeohookeanPS {
  
    NeohookeanPS() {
        double youngsModulus = 2e6;
        double poissonsRatio = 0.45;
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
    
    template<typename DataType>
    inline Eigen::Matrix33x<DataType> hessian(Eigen::MatrixBase<DataType> &stretches) {
        
        //Not implemented yet
        
    }

    double m_C, m_D;
    
};

//build specific principal stretch material
template<typename DataType, typename ShapeFunction>
using  EnergyPSNH = EnergyPrincipalStretch<DataType, ShapeFunction, NeohookeanPS>;

/* Tetrahedral finite elements */
template<typename DataType>
using FEMPSNHTet = FEMPrincipalStretchTet<DataType, EnergyPSNH>;

typedef PhysicalSystemFEM<double, FEMPSNHTet> FEMLinearTets;
//typedef PhysicalSystemFEM<double, NeohookeanTet> FEMLinearTets;

typedef World<double, std::tuple<FEMLinearTets *,PhysicalSystemParticleSingle<double> *>,
std::tuple<ForceSpringFEMParticle<double> *, ForceParticlesGravity<double> *>,
std::tuple<ConstraintFixedPoint<double> *> > MyWorld;
typedef TimeStepperEulerImplicitBFGS<double, AssemblerParallel<double, AssemblerEigenSparseMatrix<double> >,
AssemblerParallel<double, AssemblerEigenVector<double> > > MyTimeStepper;

//typedef TimeStepperEulerImplicitBFGS<double, AssemblerEigenSparseMatrix<double>,
//AssemblerEigenVector<double> > MyTimeStepper;

typedef Scene<MyWorld, MyTimeStepper> MyScene;


void preStepCallback(MyWorld &world) {
    // This is an example callback
}


//Test newton solver by solving for static equilibrium of a bendy bar
int main(int argc, char **argv) {
    std::cout<<"Test Neohookean FEM \n";
    
    //Setup Physics
    MyWorld world;
    
    //new code -- load tetgen files
    /*Eigen::MatrixXd V(4,3);
    Eigen::MatrixXi F(1,4);
    
    V << 0,0,0,
     1,0,0,
     0,1,0,
     0,0,1;
    
    F << 0, 1, 2, 3;*/
    
   /* Testing SVD derivative
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    
    readTetgen(V, F, dataDir()+"/meshesTetgen/Beam/Beam.node", dataDir()+"/meshesTetgen/Beam/Beam.ele");
    
    FEMLinearTets *test = new FEMLinearTets(V,F);
    world.addSystem(test);
    world.finalize(); //After this all we're ready to go (clean up the interface a bit later)
    
    auto q = mapStateEigen(world);
    q.setZero();
    
    std::cout << "Starting now." << std::endl;
    Eigen::VectorXi indices = minVertices(test, 0);
    Eigen::SparseMatrix<double> P = fixedPointProjectionMatrix(indices, *test,world);
    
    MyTimeStepper stepper(0.1, P);
    
    //Display
    QGuiApplication app(argc, argv);
    
    MyScene *scene = new MyScene(&world, &stepper, preStepCallback);
    GAUSSVIEW(scene);
    
    std::cout << "Done." << std::endl;
    return app.exec();*/
    
    Eigen::Vector3x<float> S;
    Eigen::Matrix33x<float> U,V;
    Eigen::Tensor3333x<float> dU, dV, dUfd, dVfd;
    Eigen::Tensor333x<float> dS, dSfd;
    
    Eigen::Matrix33x<float> A = Eigen::Matrix33x<float>::Random();
    Eigen::Matrix33x<float> Ap1, Am1;
    
    igl::svd3x3(A,U,S,V);
    
    std::cout<<"V:\n"<<V<<"\n";
    std::cout<<"S:\n"<<S.transpose()<<"\n";
    
    Eigen::dSVD(dU, dS,dV,U,S,V); 
    
    std::cout<<"A: \n"<<A<<"\n";
    
    //Finite Difference Approximations
    float dif = 1e-3;
    for(unsigned int r=0; r<3; ++r) {
        for(unsigned int s = 0; s<3; ++s) {
            Ap1 = A;
            Am1 = A;
            Ap1(r,s) += dif;
            Am1(r,s) -= dif;
    
            igl::svd3x3(Ap1,dUfd[r][s],dSfd[r][s],dVfd[r][s]);
            igl::svd3x3(Am1,U,S,V);
            dUfd[r][s] -= U;
            dUfd[r][s] /= (2.0*dif);
            
            dSfd[r][s] -= S;
            dSfd[r][s] /= (2.0*dif);
            
            dVfd[r][s] -= V;
            dVfd[r][s] /= (2.0*dif);
    
        }
    }
    
    std::cout<<"dV(0,0):\n"<<dV[1][0]<<"\n";
    std::cout<<"dVfd(0,0):\n"<<dVfd[1][0]<<"\n";
    std::cout<<"dU(0,0):\n"<<dU[1][0]<<"\n";
    std::cout<<"dUfd(0,0):\n"<<dUfd[1][0]<<"\n";
    std::cout<<"dS(0,0): "<<dS[1][0].transpose()<<"\n";
    std::cout<<"dSfd(0,0): "<<dSfd[1][0].transpose()<<"\n";

}
