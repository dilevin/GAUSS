//
//  EigenFit.h
//  Gauss
//
//  Created by Edwin Chen on 2018-05-11.
//
//

 #ifndef EigenFit_h
 #define EigenFit_h


#include <FEMIncludes.h>
#include <GaussIncludes.h>
#include <UtilitiesFEM.h>
#include <State.h>
#include <ParticleSystemIncludes.h>
#include <ConstraintFixedPoint.h>

using namespace Gauss;
using namespace FEM;
using namespace ParticleSystem; //For Force Spring

// subclass a hard-coded templated class from PhysicalSystemFEM
// this means that this EigenFit only works for NeohookeanTets
class EigenFit: public PhysicalSystemFEM<double, NeohookeanTet>{
    
public:
    // alias the hard-coded template name. Easier to read
    // the following lines read: the Physical System Implementation used here is a neo-hookean tet class
    using PhysicalSystemImpl = PhysicalSystemFEM<double, NeohookeanTet>;
    
    // use all the default function for now
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
    // the constructor will take the two mesh parameters, one coarse one fine.
    // The coarse mesh data will be passed to the parent class constructor to constructor
    // the fine mesh data will be used to initialize the members specific to the EigenFit class
     EigenFit(Eigen::MatrixXx<double> &Vc, Eigen::MatrixXi &Fc,Eigen::MatrixXx<double> &Vf, Eigen::MatrixXi &Ff, unsigned int numModes) : PhysicalSystemImpl(Vc,Fc)
     {
         // setup the fine mesh
         PhysicalSystemImpl *m_fineMeshSystem = new PhysicalSystemImpl(Vf,Ff);
         m_world.addSystem(m_fineMeshSystem);
//
         fixDisplacementMin(m_world, m_fineMeshSystem);
         m_world.finalize();
//
         auto q = mapStateEigen(m_world);
         q.setZero();
//
         // hard-coded constraint projection
         Eigen::VectorXi indices = minVertices(m_fineMeshSystem, 0);
         m_P = fixedPointProjectionMatrix(indices, *m_fineMeshSystem,m_world);

         m_numModes = numModes;
         
         m_R.setConstant(m_numModes, 0.1);
         m_I.setConstant(m_numModes, 1.0);

         
         //lambda can't capture member variable, so create a local one for lambda in ASSEMBLELIST
         // world name must match "world"?!
         World<double, std::tuple<PhysicalSystemImpl *>,
         std::tuple<ForceSpringFEMParticle<double> *, ForceParticlesGravity<double> *>,
         std::tuple<ConstraintFixedPoint<double> *> > &world = m_world;

         // assemble mass matrix in the constructor because it won't change
         
         //lambda can't capture member variable, so create a local one for lambda in ASSEMBLELIST
         AssemblerEigenSparseMatrix<double> &massMatrix = m_fineMassMatrix;
         
         //get mass matrix
         ASSEMBLEMATINIT(massMatrix, world.getNumQDotDOFs(), world.getNumQDotDOFs());
         ASSEMBLELIST(massMatrix, world.getSystemList(), getMassMatrix);
         ASSEMBLEEND(massMatrix);
         
         //constraint Projection
         (*massMatrix) = m_P*(*massMatrix)*m_P.transpose();

         
//         delete m_fineMeshSystem;
     }
    

    ~EigenFit() {
    }
    
//    void calculateEigenFitData(Eigen::MatrixXd &Y, Eigen::MatrixXd &Z, std::pair<Eigen::MatrixXx<double>, Eigen::VectorXx<double> > &m_Us){
        void calculateEigenFitData(Eigen::MatrixXd &Y, Eigen::MatrixXd &Z, std::pair<Eigen::MatrixXx<double>, Eigen::VectorXx<double> > &m_fineUs, std::pair<Eigen::MatrixXx<double>, Eigen::VectorXx<double> > &m_coarseUs){

//            std::pair<Eigen::MatrixXx<double>, Eigen::VectorXx<double> > m_Us;
        //lambda can't capture member variable, so create a local one for lambda in ASSEMBLELIST
        // world name must match "world"?!
        World<double, std::tuple<PhysicalSystemImpl *>,
        std::tuple<ForceSpringFEMParticle<double> *, ForceParticlesGravity<double> *>,
        std::tuple<ConstraintFixedPoint<double> *> > &world = m_world;
        
//        lambda can't capture member variable, so create a local one for lambda in ASSEMBLELIST
        AssemblerEigenSparseMatrix<double> &fineStiffnessMatrix = m_fineStiffnessMatrix;
        
        //get stiffness matrix
        ASSEMBLEMATINIT(fineStiffnessMatrix, world.getNumQDotDOFs(), world.getNumQDotDOFs());
        ASSEMBLELIST(fineStiffnessMatrix, world.getSystemList(), getStiffnessMatrix);
        ASSEMBLELIST(fineStiffnessMatrix, world.getForceList(), getStiffnessMatrix);
        ASSEMBLEEND(fineStiffnessMatrix);
        
        //constraint Projection
        (*fineStiffnessMatrix) = m_P*(*fineStiffnessMatrix)*m_P.transpose();

        
        //Eigendecomposition
        m_Us = generalizedEigenvalueProblem((*fineStiffnessMatrix), (*m_fineMassMatrix), m_numModes, 0.0);
        Y = (*m_fineMassMatrix)*m_Us.first*(m_R-m_I).asDiagonal();
        Z =  (m_Us.second.asDiagonal()*m_Us.first.transpose()*(*m_fineMassMatrix));
        
    }
    
    
    
protected:
    
    World<double, std::tuple<PhysicalSystemImpl *>,
    std::tuple<ForceSpringFEMParticle<double> *, ForceParticlesGravity<double> *>,
    std::tuple<ConstraintFixedPoint<double> *> > m_world;
    
    AssemblerEigenSparseMatrix<double> m_fineMassMatrix;
    AssemblerEigenSparseMatrix<double> m_fineStiffnessMatrix;

    AssemblerEigenSparseMatrix<double> m_coarseMassMatrix;
    AssemblerEigenSparseMatrix<double> m_coarseStiffnessMatrix;

    //num modes to correct
    unsigned int m_numModes;
//    unsigned int m_numToCorrect;
    
    //Ratios diagonal matrix, stored as vector
    Eigen::VectorXd m_R;
    Eigen::VectorXd m_I;
    
    //Subspace Eigenvectors and eigenvalues
    AssemblerEigenVector<double> m_fExt;
    Eigen::SparseMatrix<double> m_P;
    
    

private:
    
};

 #endif /* EigenFit_h */

