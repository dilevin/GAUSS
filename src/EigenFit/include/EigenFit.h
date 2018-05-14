//
//  EigenFit.h
//  Gauss
//
//  Created by Edwin Chen on 2018-05-11.
//
//

// #ifndef EigenFit_h
// #define EigenFit_h


// #endif /* EigenFit_h */

#include <FEMIncludes.h>
#include <GaussIncludes.h>
#include <UtilitiesFEM.h>
#include <State.h>

using namespace Gauss;
using namespace FEM;
using namespace ParticleSystem; //For Force Spring

//using FEMLinearTets = PhysicalSystemFEM<double, NeohookeanTet>;


class EigenFit: public PhysicalSystemFEM<double, NeohookeanTet>{
    
public:
    using PhysicalSystemImpl = PhysicalSystemFEM<double, NeohookeanTet>;
    using PhysicalSystemImpl::getEnergy;
    using PhysicalSystemImpl::getStrainEnergy;
    using PhysicalSystemImpl::getStrainEnergyPerElement;
    using PhysicalSystemImpl::getMassMatrix;
    using PhysicalSystemImpl::getStiffnessMatrix;
    using PhysicalSystemImpl::getForce;
    using PhysicalSystemImpl::getInternalForce;
    using PhysicalSystemImpl::getQ;
    using PhysicalSystemImpl::getQDot;
    using PhysicalSystemImpl::getPosition;
    using PhysicalSystemImpl::getDPDQ;
    using PhysicalSystemImpl::getVelocity;
    using PhysicalSystemImpl::getDVDQ;
    using PhysicalSystemImpl::getGeometry;
    
    // constructor
     EigenFit(Eigen::MatrixXx<double> &Vc, Eigen::MatrixXi &Fc,Eigen::MatrixXx<double> &Vf, Eigen::MatrixXi &Ff) : PhysicalSystemImpl(Vc,Fc)
     {
         m_fineMesh = new PhysicalSystemImpl(Vf,Ff);
         m_state = new State<double>(0,m_fineMesh->getQ().getNumScalarDOF() + m_fineMesh->getQDot().getNumScalarDOF());
         
     }
    

    // pass through constructor without fine mesh
//    EigenFit(Eigen::MatrixXx<double> &Vc, Eigen::MatrixXi &Fc) : PhysicalSystemFEM<double, NeohookeanTet>(Vc,Fc)
//    {
//    }
    

    ~EigenFit() {delete m_fineMesh; delete m_state;}
    /*
    void getStiffness(Matrix &assembler, State &state){
        // Get fine mass matrix
        auto q = //mapState(m_state)
        q = //embedding function
        
        //Get fine stiffness matrix
    }
    
    void getForce(){
        
    }
     */
    
protected:
PhysicalSystemFEM<double,NeohookeanTet> *m_fineMesh;
State<double> *m_state;
    AssemblerEigenSparseMatrix<double> m_massMatrix;
    AssemblerEigenSparseMatrix<double> m_stiffnessMatrix;

private:
    
};
