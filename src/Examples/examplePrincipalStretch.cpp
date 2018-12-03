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
    return app.exec();
    
    
}
