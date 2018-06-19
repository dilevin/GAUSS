//
//  EigenFit.h
//  Gauss
//
//  Created by Edwin Chen on 2018-05-11.
//
//
//#define EDWIN_DEBUG

#ifndef EigenFit_h
#define EigenFit_h


#include <FEMIncludes.h>
#include <GaussIncludes.h>
#include <UtilitiesFEM.h>
#include <State.h>
#include <ParticleSystemIncludes.h>
#include <ConstraintFixedPoint.h>
#include <unsupported/Eigen/SparseExtra>
#include <sys/stat.h>

using namespace Gauss;
using namespace FEM;
using namespace ParticleSystem; //For Force Spring

// subclass a hard-coded templated class from PhysicalSystemFEM
// this means that this EigenFit only works for NeohookeanHFixedTets
class EigenFit: public PhysicalSystemFEM<double, NeohookeanHFixedTet>{
    
public:
    // alias the hard-coded template name. Easier to read
    // the following lines read: the Physical System Implementation used here is a neo-hookean tet class
    using PhysicalSystemImpl = PhysicalSystemFEM<double, NeohookeanHFixedTet>;
    
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
    EigenFit(Eigen::MatrixXx<double> &Vc, Eigen::MatrixXi &Fc,Eigen::MatrixXx<double> &Vf, Eigen::MatrixXi &Ff, bool flag, double youngs, double poisson, int constraint_dir, double constraint_tol, unsigned int cswitch, unsigned int hausdorff_dist, unsigned int numModes) : PhysicalSystemImpl(Vc,Fc)
    {
        m_Vf = Vf;
        m_Ff = Ff;
        
        ratio_recalculation_flag = flag;
        constraint_switch = cswitch;
        
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
        
        
        // set up material parameters
        for(unsigned int iel=0; iel<m_fineMeshSystem->getImpl().getF().rows(); ++iel) {
            
            m_fineMeshSystem->getImpl().getElement(iel)->setParameters(youngs, poisson);
            
        }
        m_fineWorld.addSystem(m_fineMeshSystem);
        
        
        //       constraints
        Eigen::SparseMatrix<double> P;
        if (constraint_switch == 0) {
            // hard-coded constraint projection
            
            m_fineWorld.finalize();
            
            P.resize(Vf.rows()*3,Vf.rows()*3);
            P.setIdentity();
            m_P = P;
            m_numConstraints = 10;
        }
        else if (constraint_switch == 1)
        {
            // default constraint
            //            fix displacement
            fixDisplacementMin(m_fineWorld, m_fineMeshSystem, constraint_dir, constraint_tol);
            
            m_fineWorld.finalize();
            // hard-coded constraint projection
            Eigen::VectorXi indices = minVertices(m_fineMeshSystem, constraint_dir, constraint_tol);
            P = fixedPointProjectionMatrix(indices, *m_fineMeshSystem,m_fineWorld);
            m_P = P;
            m_numConstraints = indices.size();
        }
        else if (constraint_switch == 2)
        {
            
            Eigen::VectorXi movingVerts = maxVerticesTol(m_fineMeshSystem, constraint_dir, constraint_tol);//indices for moving parts
            std::vector<ConstraintFixedPoint<double> *> movingConstraints;

            for(unsigned int ii=0; ii<movingVerts.rows(); ++ii) {
                movingConstraints.push_back(new ConstraintFixedPoint<double>(&m_fineMeshSystem->getQ()[movingVerts[ii]], Eigen::Vector3d(0,0,0)));
                m_fineWorld.addConstraint(movingConstraints[ii]);
            }
            m_fineWorld.finalize(); //After this all we're ready to go (clean up the interface a bit later)
            
            // hard-coded constraint projection
            P = fixedPointProjectionMatrix(movingVerts, *m_fineMeshSystem,m_fineWorld);
            m_P = P;
            m_numConstraints = movingVerts.size();
        }
        
        
        
        //lambda can't capture member variable, so create a local one for lambda in ASSEMBLELIST
        // world name must match "world"?!
        World<double, std::tuple<PhysicalSystemImpl *>,
        std::tuple<ForceSpringFEMParticle<double> *, ForceParticlesGravity<double> *>,
        std::tuple<ConstraintFixedPoint<double> *> > &world = m_fineWorld;
        
        auto q = mapStateEigen(m_fineWorld);
        q.setZero();
        
        // the first few ratios are 1 if less than 6 constraints, because eigenvalues ratio 0/0 is not defined
        if (m_numConstraints > 6) {
            // if constraint is more than a point constaint
            m_numModes = numModes;
            
            // put random value to m_R for now
            m_R.setConstant(m_numModes, 0.1);
            ratio_calculated = false;
            m_I.setConstant(m_numModes, 1.0);
        }
        else if (m_numConstraints == 3)
        {
            // if constraint is  a point constaint
            m_numModes = numModes;
            m_numModes = m_numModes + 3;
            
            // put random value to m_R for now
            m_R.setConstant(m_numModes, 0.1);
            m_R(0) = 1.0;
            m_R(1) = 1.0;
            m_R(2) = 1.0;
            ratio_calculated = false;
            m_I.setConstant(m_numModes, 1.0);
        }
        else
        {
            // otherwise, free boundary
            m_numModes = numModes;
            m_numModes = m_numModes + 6;
            
            // put random value to m_R for now
            m_R.setConstant(m_numModes, 0.1);
            m_R(0) = 1.0;
            m_R(1) = 1.0;
            m_R(2) = 1.0;
            m_R(3) = 1.0;
            m_R(4) = 1.0;
            m_R(5) = 1.0;

            ratio_calculated = false;
            m_I.setConstant(m_numModes, 1.0);
            
        }
        
        
        
        // assemble mass matrix in the constructor because it won't change
        
        //lambda can't capture member variable, so create a local one for lambda in ASSEMBLELIST
        AssemblerEigenSparseMatrix<double> &massMatrix = m_fineMassMatrix;
        
        //get mass matrix
        ASSEMBLEMATINIT(massMatrix, world.getNumQDotDOFs(), world.getNumQDotDOFs());
        ASSEMBLELIST(massMatrix, world.getSystemList(), getMassMatrix);
        ASSEMBLEEND(massMatrix);
        
        //constraint Projection
        (*massMatrix) = m_P*(*massMatrix)*m_P.transpose();
        
        // fill in the rest state position
        restFineState = m_fineWorld.getState();
        
        
        // create a deep copy for the rest state position
        fine_pos0 = new double[world.getNumQDOFs()];
        
        Eigen::Map<Eigen::VectorXd> eigen_fine_pos0(fine_pos0,world.getNumQDOFs());
        
        Eigen::Vector3x<double> pos0;
        
        unsigned int idx;
        idx = 0;
        
        for(unsigned int vertexId=0;  vertexId < std::get<0>(m_fineWorld.getSystemList().getStorage())[0]->getGeometry().first.rows(); ++vertexId) {
            
            // because getFinePosition is in EigenFit, not another physical system Impl, so don't need getImpl()
            pos0 = this->getFinePosition(restFineState, vertexId);
            
            eigen_fine_pos0(idx) = pos0[0];
            idx++;
            eigen_fine_pos0(idx) = pos0[1];
            idx++;
            eigen_fine_pos0(idx) = pos0[2];
            idx++;
        }
        //
        
    }
    
    
    ~EigenFit() {delete fine_pos0;
    }
    
    // calculate data, TODO: the first two parameter should be const
    template<typename MatrixAssembler>
    void calculateEigenFitData(State<double> &state, MatrixAssembler &coarseMassMatrix, MatrixAssembler &coarseStiffnessMatrix,  std::pair<Eigen::MatrixXx<double>, Eigen::VectorXx<double> > &m_coarseUs, Eigen::MatrixXd &Y, Eigen::MatrixXd &Z){
        
        
        // Eigendecomposition for the coarse mesh
        m_coarseUs = generalizedEigenvalueProblem((*coarseStiffnessMatrix), (*coarseMassMatrix), m_numModes, 1e-3);
        
        if(ratio_recalculation_flag || (!ratio_calculated)){
            //            std::pair<Eigen::MatrixXx<double>, Eigen::VectorXx<double> > m_Us;
            //lambda can't capture member variable, so create a local one for lambda in ASSEMBLELIST
            // world name must match "world"?!
            World<double, std::tuple<PhysicalSystemImpl *>,
            std::tuple<ForceSpringFEMParticle<double> *, ForceParticlesGravity<double> *>,
            std::tuple<ConstraintFixedPoint<double> *> > &world = m_fineWorld;
            
            Eigen::Map<Eigen::VectorXd> fine_q = mapStateEigen<0>(m_fineWorld);
            
            //            double pd_fine_pos[world.getNumQDOFs()]; // doesn't work for MSVS
            double* pd_fine_pos = new double[world.getNumQDOFs()];
            Eigen::Map<Eigen::VectorXd> fine_pos(pd_fine_pos,world.getNumQDOFs());
            Eigen::Map<Eigen::VectorXd> eigen_fine_pos0(fine_pos0,world.getNumQDOFs());
            
            Eigen::Vector3x<double> pos;
            unsigned int idx;
            idx = 0;
            
            for(unsigned int vertexId=0;  vertexId < std::get<0>(m_fineWorld.getSystemList().getStorage())[0]->getGeometry().first.rows(); ++vertexId) {
                
                // because getFinePosition is in EigenFit, not another physical system Impl, so don't need getImpl()
                pos = this->getFinePosition(state, vertexId);
                
                fine_q(idx) = pos[0] - eigen_fine_pos0[idx];
                idx++;
                fine_q(idx) = pos[1] - eigen_fine_pos0[idx];
                idx++;
                fine_q(idx) = pos[2] - eigen_fine_pos0[idx];
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
            m_Us = generalizedEigenvalueProblem((*fineStiffnessMatrix), (*m_fineMassMatrix), m_numModes, 1e-3);
            
            for(int i = 0; i < m_numModes; ++i)
            {
                m_R(i) = m_Us.second(i)/m_coarseUs.second(i);
                if (m_numConstraints == 3)
                {
                    // if constraint is  a point constaint
                    m_R(0) = 1.0;
                    m_R(1) = 1.0;
                    m_R(2) = 1.0;
                }
                else if(m_numConstraints == 0)
                {
                    //  free boundary
                    m_R(0) = 1.0;
                    m_R(1) = 1.0;
                    m_R(2) = 1.0;
                    m_R(3) = 1.0;
                    m_R(4) = 1.0;
                    m_R(5) = 1.0;
                    
                }
#ifdef EDWIN_DEBUG
                std::cout<<m_R(i)<<std::endl;
#endif
                
            }
            ratio_calculated = true;
        }
        Y = (*coarseMassMatrix)*m_coarseUs.first*(m_R-m_I).asDiagonal();
        Z =  (m_coarseUs.second.asDiagonal()*m_coarseUs.first.transpose()*(*coarseMassMatrix));
        
    }
    
    //per vertex accessors. takes the state of the coarse mesh
    inline Eigen::Vector3x<double> getFinePosition(const State<double> &state, unsigned int vertexId) const {
        return m_Vf.row(vertexId).transpose() + m_N.block(3*vertexId, 0, 3, m_N.cols())*(*this).getImpl().getElement(m_elements[vertexId])->q(state);
    }
    
    inline World<double, std::tuple<PhysicalSystemImpl *>,
    std::tuple<ForceSpringFEMParticle<double> *, ForceParticlesGravity<double> *>,
    std::tuple<ConstraintFixedPoint<double> *> > & getFineWorld(){ return m_fineWorld;}
    
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
    
    double* fine_pos0  = NULL;
    // rest state of fine q
    
    bool ratio_recalculation_flag;
    bool ratio_calculated;
    
    unsigned int constraint_switch;
    unsigned int m_numConstraints;
private:
    
};

#endif /* EigenFit_h */

