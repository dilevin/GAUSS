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
//    using PhysicalSystemImpl::getElements;
    
    // constructor
    // the constructor will take the two mesh parameters, one coarse one fine.
    // The coarse mesh data will be passed to the parent class constructor to constructor
    // the fine mesh data will be used to initialize the members specific to the EigenFit class
     EigenFit(Eigen::MatrixXx<double> &Vc, Eigen::MatrixXi &Fc,Eigen::MatrixXx<double> &Vf, Eigen::MatrixXi &Ff) : PhysicalSystemImpl(Vc,Fc)
     {
         m_Vf = Vf;
         m_Ff = Ff;
         
         AssemblerEigenSparseMatrix<double> N;
         //element[i] is a n-vector that stores the index of the element containing the ith vertex in the embedded mesh
         getShapeFunctionMatrix(N,m_elements,Vf, (*this).getImpl());
         
         Eigen::Vector3x<double> vertex = m_Vf.row(0);
         
         unsigned int numCols = (*this).getImpl().getElements()[0]->N(vertex.data()).cols();
         unsigned int el;
         // calculate the generalized barycentric coordinates
         m_N.resize(3*m_Vf.rows(), numCols);
         for(unsigned int ii=0; ii<m_elements.size(); ++ii) {
             el = m_elements[ii];
             vertex = m_Vf.row(ii);
             
             auto Jmat = (*this).getImpl().getElements()[el]->N(vertex.data());
             
             m_N.block(3*ii, 0, 3, numCols) = Jmat;
         }

         
         // setup the fine mesh
         PhysicalSystemImpl *m_fineMeshSystem = new PhysicalSystemImpl(Vf,Ff);
         m_fineWorld.addSystem(m_fineMeshSystem);
//
         fixDisplacementMin(m_fineWorld, m_fineMeshSystem);
         m_fineWorld.finalize();
//
         auto q = mapStateEigen(m_fineWorld);
         q.setZero();
//
         // hard-coded constraint projection
         Eigen::VectorXi indices = minVertices(m_fineMeshSystem, 0);
         m_P = fixedPointProjectionMatrix(indices, *m_fineMeshSystem,m_fineWorld);

         m_numModes = 10;
         
         m_R.setConstant(m_numModes, 0.1);
         m_I.setConstant(m_numModes, 1.0);

         
         //lambda can't capture member variable, so create a local one for lambda in ASSEMBLELIST
         // world name must match "world"?!
         World<double, std::tuple<PhysicalSystemImpl *>,
         std::tuple<ForceSpringFEMParticle<double> *, ForceParticlesGravity<double> *>,
         std::tuple<ConstraintFixedPoint<double> *> > &world = m_fineWorld;

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
         
         restFineState = m_fineWorld.getState();
//
         
     }
    

    ~EigenFit() {
    }
    
//    void calculateEigenFitData(Eigen::MatrixXd &Y, Eigen::MatrixXd &Z, std::pair<Eigen::MatrixXx<double>, Eigen::VectorXx<double> > &m_Us){
    // calculate data, TODO: the first two parameter should be const
        void calculateEigenFitData(const State<double> &state, const AssemblerEigenSparseMatrix<double> &coarseMassMatrix, const AssemblerEigenSparseMatrix<double> &coarseStiffnessMatrix,  std::pair<Eigen::MatrixXx<double>, Eigen::VectorXx<double> > &m_coarseUs, Eigen::MatrixXd &Y, Eigen::MatrixXd &Z){

//            std::pair<Eigen::MatrixXx<double>, Eigen::VectorXx<double> > m_Us;
        //lambda can't capture member variable, so create a local one for lambda in ASSEMBLELIST
        // world name must match "world"?!
        World<double, std::tuple<PhysicalSystemImpl *>,
        std::tuple<ForceSpringFEMParticle<double> *, ForceParticlesGravity<double> *>,
        std::tuple<ConstraintFixedPoint<double> *> > &world = m_fineWorld;
            
            Eigen::Map<Eigen::VectorXd> fine_q = mapStateEigen<0>(m_fineWorld);

//            update the embedding q
            Eigen::Vector3x<double> pos;
            Eigen::Vector3x<double> pos0;
            unsigned int idx;
            idx = 0;
            
            for(unsigned int vertexId=0;  vertexId < std::get<0>(m_fineWorld.getSystemList().getStorage())[0]->getGeometry().first.rows(); ++vertexId) {
                
//                    vertexId = this->getImpl().getGeometry().second(elId, ii);
                    // because getFinePosition is in EigenFit, not another physical system Impl, so don't need getImpl()
                    pos = this->getFinePosition(state, vertexId);
                    pos0 = this->getFinePosition(restFineState, vertexId);
                
                    fine_q(idx) = pos[0] - pos0[0];
                    idx++;
                    fine_q(idx) = pos[1] - pos0[1];
                    idx++;
                    fine_q(idx) = pos[2] - pos0[2];
                    idx++;
                }
                
        

        
//        lambda can't capture member variable, so create a local one for lambda in ASSEMBLELIST
        AssemblerEigenSparseMatrix<double> &fineStiffnessMatrix = m_fineStiffnessMatrix;
        
        //get stiffness matrix
        ASSEMBLEMATINIT(fineStiffnessMatrix, world.getNumQDotDOFs(), world.getNumQDotDOFs());
        ASSEMBLELIST(fineStiffnessMatrix, world.getSystemList(), getStiffnessMatrix);
        ASSEMBLELIST(fineStiffnessMatrix, world.getForceList(), getStiffnessMatrix);
        ASSEMBLEEND(fineStiffnessMatrix);
        
        //constraint Projection
        (*fineStiffnessMatrix) = m_P*(*fineStiffnessMatrix)*m_P.transpose();

        //Eigendecomposition for the embedded fine mesh
        std::pair<Eigen::MatrixXx<double>, Eigen::VectorXx<double> > m_Us;
        m_Us = generalizedEigenvalueProblem((*fineStiffnessMatrix), (*m_fineMassMatrix), m_numModes, 0.0);
            
        // Eigendecomposition for the coarse mesh
        m_coarseUs = generalizedEigenvalueProblem((*coarseStiffnessMatrix), (*coarseMassMatrix), m_numModes, 0.0);
            
        for(int i = 0; i < m_numModes; ++i)
            {
                m_R(i) = m_Us.second(i)/m_coarseUs.second(i);
            }
        Y = (*coarseMassMatrix)*m_coarseUs.first*(m_R-m_I).asDiagonal();
        Z =  (m_coarseUs.second.asDiagonal()*m_coarseUs.first.transpose()*(*coarseMassMatrix));
        
    }
    
    //per vertex accessors. takes the state of the coarse mesh
    inline Eigen::Vector3x<double> getFinePosition(const State<double> &state, unsigned int vertexId) const {
        //return m_Vf.row(vertexId).transpose() + m_N.block(3*vertexId, 0, 3, m_N.cols())*mapDOFEigen(fem.getQ(), state);
        return m_Vf.row(vertexId).transpose() + m_N.block(3*vertexId, 0, 3, m_N.cols())*(*this).getImpl().getElement(m_elements[vertexId])->q(state);
    }
    
    
protected:
    
    World<double, std::tuple<PhysicalSystemImpl *>,
    std::tuple<ForceSpringFEMParticle<double> *, ForceParticlesGravity<double> *>,
    std::tuple<ConstraintFixedPoint<double> *> > m_fineWorld;
    
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
    
    Eigen::MatrixXx<double> m_Vf;
    Eigen::MatrixXi m_Ff;
    Eigen::MatrixXd m_N;
    //m_elements[i] is a n-vector that stores the index of the element containing the ith vertex in the embedded mesh
    Eigen::VectorXi m_elements;
    
    State<double> restFineState;
    // rest state of fine q
//    Eigen::Map<Eigen::VectorXd> fine_pos0;

private:
    
};

 #endif /* EigenFit_h */

