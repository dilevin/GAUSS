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
#include <igl/sortrows.h>
#include <igl/histc.h>
#include <igl/unique_rows.h>
#include <igl/hausdorff.h>
#include <fstream>
#include <SolverPardiso.h>
#include <SparseRegularInversePardiso.h>

using namespace Gauss;
using namespace FEM;
using namespace ParticleSystem; //For Force Spring
using std::cout;
using std::endl;

// subclass a hard-coded templated class from PhysicalSystemFEM
// this means that this EigenFit only works for NeohookeanHFixedTets
class EigenFit: public PhysicalSystemFEM<double, NeohookeanHFixedTet>{
    //class EigenFit: public PhysicalSystemFEM<double, NeohookeanHFixedTet>{
    
public:
    // alias the hard-coded template name. Easier to read
    // the following lines read: the Physical System Implementation used here is a neo-hookean tet class
    //    using PhysicalSystemImpl = PhysicalSystemFEM<double, NeohookeanHFixedTet>;
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
    EigenFit(Eigen::MatrixXx<double> &Vc, Eigen::MatrixXi &Fc,Eigen::MatrixXx<double> &Vf, Eigen::MatrixXi &Ff, int dynamic_switch, double youngs, double poisson, int constraintDir, double constraintTol, unsigned int cswitch, unsigned int hausdorff_dist, unsigned int numModes, std::string cmeshname, std::string fmeshname ) : PhysicalSystemImpl(Vc,Fc)
    {
        
        cout<<"haus:"<<hausdorff_dist<<endl;
        if(numModes != 0)
        {
            m_Vf = Vf;
            m_Ff = Ff;
            
            m_constraintDir = constraintDir;
            m_constraintTol = constraintTol;
            
            step_number = 0;
            
            std::cout<<m_Vf.rows()<<std::endl;
            std::cout<<m_Vf.cols()<<std::endl;
            
            ratio_recalculation_switch = dynamic_switch;
            constraint_switch = cswitch;
            
            //element[i] is a n-vector that stores the index of the element containing the ith vertex in the embedded mesh
            // *N is the upsample operator
            // (*N).transpose is downsample operator
            getShapeFunctionMatrix(N,m_elements,Vf, (*this).getImpl());
            
            Eigen::Vector3x<double> vertex = m_Vf.row(0);
            
            // col of the shape func, I think. 12 for Tet, 24 for hex
            unsigned int numCols = (*this).getImpl().getElements()[0]->N(vertex.data()).cols();
            unsigned int el;
            
            // set the flag
            haus = hausdorff_dist;
            
            // setup the fine mesh
            PhysicalSystemImpl *m_fineMeshSystem = new PhysicalSystemImpl(Vf,Ff);
            
            
            // set up material parameters
            this->youngs = youngs;
            this->poisson = poisson;
            for(unsigned int iel=0; iel<m_fineMeshSystem->getImpl().getF().rows(); ++iel) {
                
                m_fineMeshSystem->getImpl().getElement(iel)->setParameters(youngs, poisson);
                
            }
            m_fineWorld.addSystem(m_fineMeshSystem);
            
            //            Eigen::Map<Eigen::VectorXd> fine_q = mapStateEigen<0>(m_fineWorld);
            //            fine_q_transfered = new double[fine_q.rows()];
            
            
            //            Eigen::Map<Eigen::VectorXd> eigen_fineq_q_transfered(fine_q_transfered,fine_q.rows());
            //            eigen_fineq_q_transfered.setZero();
            
            //            fine_q_transfered = mapStateEigen<0>(m_fineWorld);
            
            //       constraints
            Eigen::SparseMatrix<double> fineP;
            Eigen::SparseMatrix<double> coarseP;
            if (constraint_switch == 0) {
                // hard-coded constraint projection
                
                m_fineWorld.finalize();
                
                fineP.resize(Vf.rows()*3,Vf.rows()*3);
                fineP.setIdentity();
                m_fineP = fineP;
                coarseP.resize(Vc.rows()*3,Vc.rows()*3);
                coarseP.setIdentity();
                m_coarseP = coarseP;
                
                m_numConstraints = 0;
            }
            else if (constraint_switch == 1)
            {
                // default constraint
                //            fix displacement
                fixDisplacementMin(m_fineWorld, m_fineMeshSystem, constraintDir, constraintTol);
                
                m_fineWorld.finalize();
                // hard-coded constraint projection
                Eigen::VectorXi fineFixedVerts = minVertices(m_fineMeshSystem, constraintDir, constraintTol);
                m_fineFixedVerts = fineFixedVerts;
                fineP = fixedPointProjectionMatrix(fineFixedVerts, *m_fineMeshSystem,m_fineWorld);
                m_fineP = fineP;
                
                Eigen::VectorXi coarse_constrait_indices = minVertices(this, constraintDir, constraintTol);
                coarseP = fixedPointProjectionMatrixCoarse(coarse_constrait_indices);
                m_coarseP = coarseP;
                
                // only need to record one because only need to know if it's 0, 3, or 6. either fine or coarse would work
                m_numConstraints = fineFixedVerts.size();
                
            }
            else if (constraint_switch == 2)
            {
                
                //                Eigen::VectorXi fineMovingVerts = minVertices(m_fineMeshSystem, constraintDir, constraintTol);//indices for moving parts
                Eigen::VectorXi fineMovingVerts;
                
                // read constraints
                
                Eigen::loadMarketVector(fineMovingVerts,  "def_init/" + fmeshname + "_fixed_min_verts.mtx");
                m_fineMovingVerts = fineMovingVerts;
                
                std::vector<ConstraintFixedPoint<double> *> fineMovingConstraints;
                
                for(unsigned int ii=0; ii<fineMovingVerts.rows(); ++ii) {
                    fineMovingConstraints.push_back(new ConstraintFixedPoint<double>(&m_fineMeshSystem->getQ()[fineMovingVerts[ii]], Eigen::Vector3d(0,0,0)));
                    m_fineWorld.addConstraint(fineMovingConstraints[ii]);
                }
                m_fineWorld.finalize(); //After this all we're ready to go (clean up the interface a bit later)
                
                // hard-coded constraint projection
                fineP = fixedPointProjectionMatrix(fineMovingVerts, *m_fineMeshSystem,m_fineWorld);
                m_fineP = fineP;
                // only need to record one because only need to know if it's 0, 3, or 6. either fine or coarse is fine
                m_numConstraints = fineMovingVerts.size();
                
                //                Eigen::VectorXi coarseMovingVerts = minVertices(this, constraintDir, constraintTol);
                Eigen::VectorXi coarseMovingVerts;
                Eigen::loadMarketVector(coarseMovingVerts,  "def_init/" + cmeshname + "_fixed_min_verts.mtx");
                
                coarseP = fixedPointProjectionMatrixCoarse(coarseMovingVerts);
                m_coarseP = coarseP;
                
            }
            else if (constraint_switch == 3)
            {
                std::vector<ConstraintFixedPoint<double> *> fixedConstraints;
                Eigen::VectorXi fineFixedVerts;
                
                // read constraints
                
                Eigen::loadMarketVector(fineFixedVerts,  "def_init/" + fmeshname + "_fixed_min_verts.mtx");
                m_fineFixedVerts = fineFixedVerts;
                //
                for(unsigned int ii=0; ii<fineFixedVerts.rows(); ++ii) {
                    fixedConstraints.push_back(new ConstraintFixedPoint<double>(&m_fineMeshSystem->getQ()[fineFixedVerts[ii]], Eigen::Vector3d(0,0,0)));
                    m_fineWorld.addConstraint(fixedConstraints[ii]);
                }
                
                m_fineWorld.finalize(); //After this all we're ready to go (clean up the interface a bit later)
                
                fineP = fixedPointProjectionMatrix(fineFixedVerts, *m_fineMeshSystem,m_fineWorld);
                m_fineP = fineP;
                
                Eigen::VectorXi coarse_constrait_indices;
                Eigen::loadMarketVector(coarse_constrait_indices,  "def_init/" + cmeshname + "_fixed_min_verts.mtx");
                
                coarseP = fixedPointProjectionMatrixCoarse(coarse_constrait_indices);
                m_coarseP = coarseP;
                
                m_numConstraints = fineFixedVerts.rows();
                
                
            }
            else if (constraint_switch == 4 || constraint_switch == 5 || constraint_switch == 6 || constraint_switch == 7 || constraint_switch == 8)
            {
                
                Eigen::VectorXi fineMovingVerts = minVertices(m_fineMeshSystem, constraintDir, constraintTol);//indices for moving parts
                m_fineMovingVerts = fineMovingVerts;
                
                std::vector<ConstraintFixedPoint<double> *> fineMovingConstraints;
                
                for(unsigned int ii=0; ii<fineMovingVerts.rows(); ++ii) {
                    fineMovingConstraints.push_back(new ConstraintFixedPoint<double>(&m_fineMeshSystem->getQ()[fineMovingVerts[ii]], Eigen::Vector3d(0,0,0)));
                    m_fineWorld.addConstraint(fineMovingConstraints[ii]);
                }
                m_fineWorld.finalize(); //After this all we're ready to go (clean up the interface a bit later)
                
                // hard-coded constraint projection
                fineP = fixedPointProjectionMatrix(fineMovingVerts, *m_fineMeshSystem,m_fineWorld);
                m_fineP = fineP;
                // only need to record one because only need to know if it's 0, 3, or 6. either fine or coarse is fine
                m_numConstraints = fineMovingVerts.size();
                
                Eigen::VectorXi coarseMovingVerts = minVertices(this, constraintDir, constraintTol);
                m_coarseMovingVerts = coarseMovingVerts;
                
                coarseP = fixedPointProjectionMatrixCoarse(coarseMovingVerts);
                m_coarseP = coarseP;
                
            }
            
            
            //lambda can't capture member variable, so create a local one for lambda in ASSEMBLELIST
            // world name must match "world"?!
            World<double, std::tuple<PhysicalSystemImpl *>,
            std::tuple<ForceSpringFEMParticle<double> *, ForceParticlesGravity<double> *>,
            std::tuple<ConstraintFixedPoint<double> *> > &world = m_fineWorld;
            
            auto q = mapStateEigen<0>(m_fineWorld);
            q.setZero();
            
            auto v = mapStateEigen<1>(m_fineWorld);
            v.setZero();
            
            // the first few ratios are 1 if less than 6 constraints, because eigenvalues ratio 0/0 is not defined
            if (m_numConstraints > 6) {
                // if constraint is more than a point constaint
                m_numModes = numModes;
                
                // put random value to m_R for now
                m_R.setConstant(m_numModes, 1.0);
                ratio_calculated = false;
                m_I.setConstant(m_numModes, 1.0);
            }
            else if (m_numConstraints == 3)
            {
                // if constraint is  a point constaint
                m_numModes = numModes;
                m_numModes = m_numModes + 3;
                
                // put random value to m_R for now
                m_R.setConstant(m_numModes, 1.0);
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
                m_R.setConstant(m_numModes, 1.0);
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
            (*massMatrix) = m_fineP*(*massMatrix)*m_fineP.transpose();
            
            // fill in the rest state position
            restFineState = m_fineWorld.getState();
            
            
            // create a deep copy for the rest state position
            fine_pos0 = new double[world.getNumQDOFs()];
            Eigen::Map<Eigen::VectorXd> eigen_fine_pos0(fine_pos0,world.getNumQDOFs());
            
            Eigen::Vector3x<double> pos0;
            
            unsigned int idx;
            idx = 0;
            
            for(unsigned int vertexId=0;  vertexId < std::get<0>(m_fineWorld.getSystemList().getStorage())[0]->getGeometry().first.rows(); ++vertexId) {
                
                pos0 = m_Vf.row(vertexId);
                
                eigen_fine_pos0(idx) = pos0[0];
                idx++;
                eigen_fine_pos0(idx) = pos0[1];
                idx++;
                eigen_fine_pos0(idx) = pos0[2];
                idx++;
            }
            
        }
    }
    
    
    ~EigenFit() {delete fine_pos0;
    }
    
    void calculateFineMesh(){
        World<double, std::tuple<PhysicalSystemImpl *>,
        std::tuple<ForceSpringFEMParticle<double> *, ForceParticlesGravity<double> *>,
        std::tuple<ConstraintFixedPoint<double> *> > &world = m_fineWorld;
        
        //        Eigen::Map<Eigen::VectorXd> fine_q = mapStateEigen<0>(m_fineWorld);
        
        //            double pd_fine_pos[world.getNumQDOFs()]; // doesn't work for MSVS
        //        Eigen::Map<Eigen::VectorXd> eigen_fine_pos0(fine_pos0,world.getNumQDOFs());
        
        //        Eigen::VectorXx<double> posFull;
        //        posFull = this->getFinePositionFull(q);
        //
        //        fine_q = posFull - eigen_fine_pos0;
        //        fine_q.setZero();
        
        //        lambda can't capture member variable, so create a local one for lambda in ASSEMBLELIST
        AssemblerEigenSparseMatrix<double> &fineStiffnessMatrix = m_fineStiffnessMatrix;
        
        //            std::cout<<
        
        //get stiffness matrix
        ASSEMBLEMATINIT(fineStiffnessMatrix, world.getNumQDotDOFs(), world.getNumQDotDOFs());
        ASSEMBLELIST(fineStiffnessMatrix, world.getSystemList(), getStiffnessMatrix);
        ASSEMBLELIST(fineStiffnessMatrix, world.getForceList(), getStiffnessMatrix);
        ASSEMBLEEND(fineStiffnessMatrix);
        //        Eigen::saveMarket(fine_q, "fineq0.dat");
        
        // for edwin debug
        AssemblerEigenVector<double> &fineforceVector = m_fineforceVector;
        AssemblerEigenVector<double> &finefExt = m_finefExt;
        
        //Need to filter internal forces seperately for this applicat
        ASSEMBLEVECINIT(fineforceVector, world.getNumQDotDOFs());
        ASSEMBLELIST(fineforceVector, world.getSystemList(), getImpl().getInternalForce);
        ASSEMBLEEND(fineforceVector);
        
        ASSEMBLEVECINIT(finefExt, world.getNumQDotDOFs());
        ASSEMBLELIST(finefExt, world.getSystemList(), getImpl().getBodyForce);
        ASSEMBLEEND(finefExt);
        
        // add external force
        (*fineforceVector) = m_fineP * (*fineforceVector);
        
        
        (*fineforceVector) = (*fineforceVector) + m_fineP *(*finefExt);
        
        //constraint Projection
        (*fineStiffnessMatrix) = m_fineP*(*fineStiffnessMatrix)*m_fineP.transpose();
        
        //Eigendecomposition for the embedded fine mesh
        std::pair<Eigen::MatrixXx<double>, Eigen::VectorXx<double> > m_Us;
        m_Us = generalizedEigenvalueProblemNegative((*fineStiffnessMatrix), (*m_fineMassMatrix), m_numModes, 0.00);
        
        fineEigMassProj = m_Us;
        fineEig = m_Us;
        fineEigMassProj.first = (*m_fineMassMatrix)*fineEigMassProj.first;
    }
    
    // calculate data, TODO: the first two parameter should be const
    template<typename MatrixAssembler>
    //    void calculateEigenFitData(State<double> &state, MatrixAssembler &coarseMassMatrix, MatrixAssembler &coarseStiffnessMatrix,  std::pair<Eigen::MatrixXx<double>, Eigen::VectorXx<double> > &m_coarseUs, Eigen::MatrixXd &Y, Eigen::MatrixXd &Z){
    bool calculateEigenFitData(const Eigen::VectorXx<double> &q, MatrixAssembler &coarseMassMatrix, MatrixAssembler &coarseStiffnessMatrix,  std::pair<Eigen::MatrixXx<double>, Eigen::VectorXx<double> > &m_coarseUs, Eigen::MatrixXd &Y, Eigen::MatrixXd &Z){
        
        // matrices passed in already eliminated the constraints
        // Eigendecomposition for the coarse mesh
        m_coarseUs = generalizedEigenvalueProblemNegative((*coarseStiffnessMatrix), (*coarseMassMatrix), m_numModes , 0.00);
        
        coarseEigMassProj = m_coarseUs;
        coarseEig = m_coarseUs;
        coarseEigMassProj.first = (*coarseMassMatrix)*coarseEigMassProj.first;
        
        
        std::cout<<"dynamic switch: "<<ratio_recalculation_switch<<std::endl;
        if((!ratio_calculated))
        {
            cout<<"recalculating ratio\n";
            if( ratio_recalculation_switch == 1 || ratio_recalculation_switch == 0)
            {
                //            Hausdorff distance check
                unsigned int mode = 0;
                unsigned int idx = 0;
                Eigen::VectorXd coarse_eig_def;
                for (mode = 0; mode < m_numModes; ++mode) {
                    coarse_eig_def = (m_coarseP.transpose()*m_coarseUs.first.col(mode)).transpose();
                    //        //
                    idx = 0;
                    //                    // getGeometry().first is V
                    Eigen::MatrixXd coarse_V_disp_p = this->getImpl().getV();
                    Eigen::MatrixXd coarse_V_disp_n = this->getImpl().getV();
                    for(unsigned int vertexId=0;  vertexId < this->getImpl().getV().rows(); ++vertexId) {
                        coarse_V_disp_p(vertexId,0) += (1*coarse_eig_def(idx));
                        coarse_V_disp_n(vertexId,0) -= (1*coarse_eig_def(idx));
                        idx++;
                        coarse_V_disp_p(vertexId,1) += (1*coarse_eig_def(idx));
                        coarse_V_disp_n(vertexId,1) -= (1*coarse_eig_def(idx));
                        idx++;
                        coarse_V_disp_p(vertexId,2) += (1*coarse_eig_def(idx));
                        coarse_V_disp_n(vertexId,2) -= (1*coarse_eig_def(idx));
                        idx++;
                    }
                    Eigen::MatrixXi coarse_F = surftri(this->getImpl().getV(), this->getImpl().getF());
                    igl::writeOBJ("coarse_mesh_eigen_mode_p" + std::to_string(mode) + ".obj" ,coarse_V_disp_p, coarse_F);
                    igl::writeOBJ("coarse_mesh_eigen_mode_n" + std::to_string(mode) + ".obj",coarse_V_disp_n, coarse_F);
                }
                
                //            std::pair<Eigen::MatrixXx<double>, Eigen::VectorXx<double> > m_Us;
                //lambda can't capture member variable, so create a local one for lambda in ASSEMBLELIST
                // world name must match "world"?!
                World<double, std::tuple<PhysicalSystemImpl *>,
                std::tuple<ForceSpringFEMParticle<double> *, ForceParticlesGravity<double> *>,
                std::tuple<ConstraintFixedPoint<double> *> > &world = m_fineWorld;
                
                Eigen::Map<Eigen::VectorXd> fine_q = mapStateEigen<0>(m_fineWorld);
                
                //            double pd_fine_pos[world.getNumQDOFs()]; // doesn't work for MSVS
                Eigen::Map<Eigen::VectorXd> eigen_fine_pos0(fine_pos0,world.getNumQDOFs());
                
                Eigen::VectorXx<double> posFull;
                posFull = this->getFinePositionFull(q);
                //
                fine_q = posFull - eigen_fine_pos0;
                
                //        lambda can't capture member variable, so create a local one for lambda in ASSEMBLELIST
                AssemblerEigenSparseMatrix<double> &fineStiffnessMatrix = m_fineStiffnessMatrix;
                
                //            std::cout<<
                
                //get stiffness matrix
                ASSEMBLEMATINIT(fineStiffnessMatrix, world.getNumQDotDOFs(), world.getNumQDotDOFs());
                ASSEMBLELIST(fineStiffnessMatrix, world.getSystemList(), getStiffnessMatrix);
                ASSEMBLELIST(fineStiffnessMatrix, world.getForceList(), getStiffnessMatrix);
                ASSEMBLEEND(fineStiffnessMatrix);
                
                // for edwin debug
                AssemblerEigenVector<double> &fineforceVector = m_fineforceVector;
                AssemblerEigenVector<double> &finefExt = m_finefExt;
                
                //Need to filter internal forces seperately for this applicat
                ASSEMBLEVECINIT(fineforceVector, world.getNumQDotDOFs());
                ASSEMBLELIST(fineforceVector, world.getSystemList(), getImpl().getInternalForce);
                ASSEMBLEEND(fineforceVector);
                
                ASSEMBLEVECINIT(finefExt, world.getNumQDotDOFs());
                ASSEMBLELIST(finefExt, world.getSystemList(), getImpl().getBodyForce);
                ASSEMBLEEND(finefExt);
                
                // add external force
                (*fineforceVector) = m_fineP * (*fineforceVector);
                
                (*fineforceVector) = (*fineforceVector) + m_fineP *(*finefExt);
                
                //constraint Projection
                (*fineStiffnessMatrix) = m_fineP*(*fineStiffnessMatrix)*m_fineP.transpose();
                
                
                //Eigendecomposition for the embedded fine mesh
                std::pair<Eigen::MatrixXx<double>, Eigen::VectorXx<double> > m_Us;
                m_Us = generalizedEigenvalueProblemNegative(((*fineStiffnessMatrix)), (*m_fineMassMatrix), m_numModes, 0.00);
                
                fineEigMassProj = m_Us;
                fineEig = m_Us;
                fineEigMassProj.first = (*m_fineMassMatrix)*fineEigMassProj.first;
                
                //            Hausdorff distance check
                //
                mode = 0;
                Eigen::VectorXd fine_eig_def;
                for (mode = 0; mode < m_numModes; ++mode) {
                    fine_eig_def = (m_fineP.transpose()*m_Us.first.col(mode)).transpose();
                    idx = 0;
                    // getGeometry().first is V
                    Eigen::MatrixXd fine_V_disp = std::get<0>(m_fineWorld.getSystemList().getStorage())[0]->getGeometry().first;
                    for(unsigned int vertexId=0;  vertexId < std::get<0>(m_fineWorld.getSystemList().getStorage())[0]->getGeometry().first.rows(); ++vertexId) {
                        
                        // because getFinePosition is in EigenFit, not another physical system Impl, so don't need getImpl()
                        fine_V_disp(vertexId,0) += (1*fine_eig_def(idx));
                        idx++;
                        fine_V_disp(vertexId,1) += (1*fine_eig_def(idx));
                        idx++;
                        fine_V_disp(vertexId,2) += (1*fine_eig_def(idx));
                        idx++;
                    }
                    
                    Eigen::MatrixXi fine_F = surftri(std::get<0>(m_fineWorld.getSystemList().getStorage())[0]->getGeometry().first, std::get<0>(m_fineWorld.getSystemList().getStorage())[0]->getGeometry().second);
                    igl::writeOBJ("fine_mesh_eigen_mode" + std::to_string(mode) + ".obj",fine_V_disp,fine_F);
                    
                    Eigen::MatrixXd coarse_V_disp_p;
                    Eigen::MatrixXd coarse_V_disp_n;
                    Eigen::MatrixXi coarse_F;
                    igl::readOBJ("coarse_mesh_eigen_mode_p" + std::to_string(mode) + ".obj",coarse_V_disp_p, coarse_F);
                    igl::readOBJ("coarse_mesh_eigen_mode_n" + std::to_string(mode) + ".obj",coarse_V_disp_n, coarse_F);
                    
                    double dist_p, dist_n, dist_scaled;
                    Eigen::MatrixXd coarse_V_disp = this->getImpl().getV();
                    
                    Eigen::Array3d xyz_scales(coarse_V_disp.col(0).maxCoeff() - coarse_V_disp.col(0).minCoeff(), coarse_V_disp.col(1).maxCoeff() - coarse_V_disp.col(1).minCoeff(),coarse_V_disp.col(2).maxCoeff() - coarse_V_disp.col(2).minCoeff());
                    double max_scale = xyz_scales.abs().maxCoeff();
                    igl::hausdorff(fine_V_disp, fine_F, coarse_V_disp_p, coarse_F, dist_p);
                    igl::hausdorff(fine_V_disp, fine_F, coarse_V_disp_n, coarse_F, dist_n);
                    if(dist_p < dist_n) dist_scaled = dist_p/max_scale;
                    else dist_scaled = dist_n/max_scale;
                    
                    // fail on hausdorff distance check
                    if( dist_scaled > 0.4 )
                    {
                        std::cout<<"dist_scaled: "<< dist_scaled<<std::endl;
                        if(haus)
                        {
                            // if required to check hausdorff and failed, return 1
                            return true;
                        }
                    }
                }
                
                //                std::cout<<"m_numConstraints "<< m_numConstraints << std::endl;
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
                    //#ifdef EDWIN_DEBUG
                    std::cout<<m_R(i)<<std::endl;
                    //#endif
                    
                }
                //                Eigen::saveMarketVector(m_Us.second, cmeshnameActual + "FullState" + std::to_string(file_ind) + ".mtx");
                ratio_calculated = true;
                
            }
            else if (ratio_recalculation_switch == 2)
            {
                // local variables for reset
                World<double, std::tuple<PhysicalSystemImpl *>,
                std::tuple<ForceSpringFEMParticle<double> *, ForceParticlesGravity<double> *>,
                std::tuple<ConstraintFixedPoint<double> *> > m_fineWorld_reset;
                
                AssemblerEigenSparseMatrix<double> m_fineStiffnessMatrix_reset;
                AssemblerEigenSparseMatrix<double> m_fineMassMatrix_reset;
                
                // get the fine mesh position to reset the mesh
                World<double, std::tuple<PhysicalSystemImpl *>,
                std::tuple<ForceSpringFEMParticle<double> *, ForceParticlesGravity<double> *>,
                std::tuple<ConstraintFixedPoint<double> *> > &world = m_fineWorld;
                
                Eigen::Map<Eigen::VectorXd> fine_q = mapStateEigen<0>(m_fineWorld);
                
                //            double pd_fine_pos[world.getNumQDOFs()]; // doesn't work for MSVS
                Eigen::Map<Eigen::VectorXd> eigen_fine_pos0(fine_pos0,world.getNumQDOFs());
                
                Eigen::VectorXx<double> posFull;
                posFull = this->getFinePositionFull(q);
                //
                fine_q = posFull - eigen_fine_pos0;
                
                //                cout<<"original q size: "<<fine_q.rows()<<endl;
                unsigned int idx = 0;
                
                Eigen::MatrixXd fine_V_disp = std::get<0>(m_fineWorld.getSystemList().getStorage())[0]->getGeometry().first;
                
                for(unsigned int vertexId=0;  vertexId < std::get<0>(m_fineWorld.getSystemList().getStorage())[0]->getGeometry().first.rows(); ++vertexId) {
                    
                    // because getFinePosition is in EigenFit, not another physical system Impl, so don't need getImpl()
                    fine_V_disp(vertexId,0) += (1*fine_q(idx));
                    idx++;
                    fine_V_disp(vertexId,1) += (1*fine_q(idx));
                    idx++;
                    fine_V_disp(vertexId,2) += (1*fine_q(idx));
                    idx++;
                }
                
                //                cout<<"reset Vf size: "<<fine_V_disp.rows()<<" "<<fine_V_disp.cols()<<endl;
                // setup the fine mesh
                
                
                PhysicalSystemImpl *m_fineMeshSystem_reset = new PhysicalSystemImpl(fine_V_disp,m_Ff);
                
                Eigen::MatrixXi fine_F = surftri(fine_V_disp, std::get<0>(m_fineWorld.getSystemList().getStorage())[0]->getGeometry().second);
                
                igl::writeOBJ("fine_mesh_eigen_rest_reset.obj",fine_V_disp,fine_F);
                igl::writeOBJ("fine_mesh_eigen_rest_original.obj",std::get<0>(m_fineWorld.getSystemList().getStorage())[0]->getGeometry().first,fine_F);
                //
                //
                //                // set up material parameters
                //                for(unsigned int iel=0; iel<m_fineMeshSystem_reset->getImpl().getF().rows(); ++iel) {
                //
                //                    m_fineMeshSystem_reset->getImpl().getElement(iel)->setParameters(youngs, poisson);
                //
                //                }
                //
                m_fineWorld_reset.addSystem(m_fineMeshSystem_reset);
                
                // do the same for coarse mesh
                World<double, std::tuple<PhysicalSystemImpl *>,
                std::tuple<ForceSpringFEMParticle<double> *, ForceParticlesGravity<double> *>,
                std::tuple<ConstraintFixedPoint<double> *> > m_coarseWorld_reset;
                
                AssemblerEigenSparseMatrix<double> m_coarseStiffnessMatrix_reset;
                AssemblerEigenSparseMatrix<double> m_coarseMassMatrix_reset;
                
                
                //                cout<<"original q size: "<<fine_q.rows()<<endl;
                idx = 0;
                
                Eigen::MatrixXd V_disp = this->getGeometry().first;
                
                for(unsigned int vertexId=0;  vertexId < this->getGeometry().first.rows(); ++vertexId) {
                    
                    V_disp(vertexId,0) += (1*q(idx));
                    idx++;
                    V_disp(vertexId,1) += (1*q(idx));
                    idx++;
                    V_disp(vertexId,2) += (1*q(idx));
                    idx++;
                }
                
                //                cout<<"reset Vc size: "<<V_disp.rows()<<" "<<V_disp.cols()<<endl;
                // setup the fine mesh
                PhysicalSystemImpl *m_coarseMeshSystem_reset = new PhysicalSystemImpl(V_disp,this->getGeometry().second);
                
                Eigen::MatrixXi coarse_F = surftri(V_disp, this->getGeometry().second);
                
                igl::writeOBJ("coarse_mesh_eigen_rest_reset.obj",V_disp,coarse_F);
                igl::writeOBJ("coarse_mesh_eigen_rest_original.obj",this->getGeometry().first,coarse_F);
                
                //
                //
                //                // set up material parameters
                //                for(unsigned int iel=0; iel<m_coarseMeshSystem_reset->getImpl().getF().rows(); ++iel) {
                //
                //                    m_coarseMeshSystem_reset->getImpl().getElement(iel)->setParameters(youngs, poisson);
                //
                //                }
                
                m_coarseWorld_reset.addSystem(m_coarseMeshSystem_reset);
                
                
                // reset the constraints and finalize the new world
                //       constraints
                if (constraint_switch == 0) {
                    // hard-coded constraint projection
                    
                    m_fineWorld_reset.finalize();
                    m_coarseWorld_reset.finalize();
                    
                }
                else if (constraint_switch == 4 || constraint_switch == 5 || constraint_switch == 6 || constraint_switch == 7 || constraint_switch == 8)
                {
                    
                    std::vector<ConstraintFixedPoint<double> *> fineMovingConstraints;
                    
                    
                    for(unsigned int ii=0; ii<m_fineMovingVerts.rows(); ++ii) {
                        fineMovingConstraints.push_back(new ConstraintFixedPoint<double>(&m_fineMeshSystem_reset->getQ()[m_fineMovingVerts[ii]], Eigen::Vector3d(0,0,0)));
                        m_fineWorld_reset.addConstraint(fineMovingConstraints[ii]);
                    }
                    m_fineWorld_reset.finalize(); //After this all we're ready to go (clean up the interface a bit later)
                    
                    std::vector<ConstraintFixedPoint<double> *> coarseMovingConstraints;
                    
                    for(unsigned int ii=0; ii<m_coarseMovingVerts.rows(); ++ii) {
                        coarseMovingConstraints.push_back(new ConstraintFixedPoint<double>(&m_coarseMeshSystem_reset->getQ()[m_coarseMovingVerts[ii]], Eigen::Vector3d(0,0,0)));
                        m_coarseWorld_reset.addConstraint(coarseMovingConstraints[ii]);
                    }
                    m_coarseWorld_reset.finalize(); //After this all we're ready to go (clean up the interface a bit later)
                    
                }
                else{
                    cout<<"constraint switch "<<constraint_switch<<" not supported for dynamic alg"<<endl;
                    exit(1);
                }
                
                // construct matrices
                //lambda can't capture member variable, so create a local one for lambda in ASSEMBLELIST
                World<double, std::tuple<PhysicalSystemImpl *>,
                std::tuple<ForceSpringFEMParticle<double> *, ForceParticlesGravity<double> *>,
                std::tuple<ConstraintFixedPoint<double> *> > &world_reset = m_fineWorld_reset;
                
                auto fine_q_reset = mapStateEigen<0>(m_fineWorld_reset);
                fine_q_reset.setZero();
                
                auto fine_v_reset = mapStateEigen<1>(m_fineWorld_reset);
                fine_v_reset.setZero();
                
                //                cout<<"reset q size: "<<fine_q_reset.rows()<<endl;
                // assemble the matrices for fine mesh
                AssemblerEigenSparseMatrix<double> &massMatrix_reset = m_fineMassMatrix_reset;
                
                //get mass matrix
                ASSEMBLEMATINIT(massMatrix_reset, world_reset.getNumQDotDOFs(), world_reset.getNumQDotDOFs());
                ASSEMBLELIST(massMatrix_reset, world_reset.getSystemList(), getMassMatrix);
                ASSEMBLEEND(massMatrix_reset);
                
                //constraint Projection
                //                cout<<"mass matrix reset size: "<<(*massMatrix_reset).rows()<<" "<<(*massMatrix_reset).cols()<<endl;
                (*massMatrix_reset) = m_fineP*(*massMatrix_reset)*m_fineP.transpose();
                
                //        lambda can't capture member variable, so create a local one for lambda in ASSEMBLELIST
                AssemblerEigenSparseMatrix<double> &fineStiffnessMatrix_reset = m_fineStiffnessMatrix_reset;
                
                //get stiffness matrix
                ASSEMBLEMATINIT(fineStiffnessMatrix_reset, world_reset.getNumQDotDOFs(), world_reset.getNumQDotDOFs());
                ASSEMBLELIST(fineStiffnessMatrix_reset, world_reset.getSystemList(), getStiffnessMatrix);
                ASSEMBLELIST(fineStiffnessMatrix_reset, world_reset.getForceList(), getStiffnessMatrix);
                ASSEMBLEEND(fineStiffnessMatrix_reset);
                
                
                //constraint Projection
                (*fineStiffnessMatrix_reset) = m_fineP*(*fineStiffnessMatrix_reset)*m_fineP.transpose();
                
                //Eigendecomposition for the embedded fine mesh
                std::pair<Eigen::MatrixXx<double>, Eigen::VectorXx<double> > m_fineUs_reset;
                m_fineUs_reset = generalizedEigenvalueProblemNegative(((*fineStiffnessMatrix_reset)), (*m_fineMassMatrix_reset), m_numModes, 0.00);
                //
//                unsigned int mode = 0;
//                Eigen::VectorXd fine_eig_def;
//                for (mode = 0; mode < m_numModes; ++mode)
//                {
                    //                    fine_eig_def = (m_fineP.transpose()*m_fineUs_reset.first.col(mode)).transpose();
                    //                    idx = 0;
                    //                    // getGeometry().first is V
                    //                    Eigen::MatrixXd fine_V_disp = std::get<0>(m_fineWorld_reset.getSystemList().getStorage())[0]->getGeometry().first;
                    //                    for(unsigned int vertexId=0;  vertexId < std::get<0>(m_fineWorld_reset.getSystemList().getStorage())[0]->getGeometry().first.rows(); ++vertexId) {
                    //
                    //                        fine_V_disp(vertexId,0) += (1*fine_eig_def(idx));
                    //                        idx++;
                    //                        fine_V_disp(vertexId,1) += (1*fine_eig_def(idx));
                    //                        idx++;
                    //                        fine_V_disp(vertexId,2) += (1*fine_eig_def(idx));
                    //                        idx++;
                    //                    }
                    
                    //                    Eigen::MatrixXi fine_F = surftri(std::get<0>(m_fineWorld_reset.getSystemList().getStorage())[0]->getGeometry().first, std::get<0>(m_fineWorld_reset.getSystemList().getStorage())[0]->getGeometry().second);
                    //                    igl::writeOBJ("fine_mesh_eigen_mode_reset" + std::to_string(mode) + ".obj",fine_V_disp,fine_F);
//                }
                    //do the same for coarse mesh
                    World<double, std::tuple<PhysicalSystemImpl *>,
                    std::tuple<ForceSpringFEMParticle<double> *, ForceParticlesGravity<double> *>,
                    std::tuple<ConstraintFixedPoint<double> *> > &coarse_world_reset = m_coarseWorld_reset;
                    
                    auto coarse_q_reset = mapStateEigen<0>(m_coarseWorld_reset);
                    coarse_q_reset.setZero();
                    
                    auto coarse_v_reset = mapStateEigen<1>(m_coarseWorld_reset);
                    coarse_v_reset.setZero();
                    
                    //                    cout<<"reset q size: "<<coarse_q_reset.rows()<<endl;
                    // assemble mass matrix in the constructor because it won't change
                    
                    //lambda can't capture member variable, so create a local one for lambda in ASSEMBLELIST
                    AssemblerEigenSparseMatrix<double> &coarseMassMatrix_reset = m_coarseMassMatrix_reset;
                    
                    //get mass matrix
                    ASSEMBLEMATINIT(coarseMassMatrix_reset, coarse_world_reset.getNumQDotDOFs(), coarse_world_reset.getNumQDotDOFs());
                    ASSEMBLELIST(coarseMassMatrix_reset, coarse_world_reset.getSystemList(), getMassMatrix);
                    ASSEMBLEEND(coarseMassMatrix_reset);
                    
                    //constraint Projection
                    //                    cout<<"mass matrix reset size: "<<(*coarseMassMatrix_reset).rows()<<" "<<(*coarseMassMatrix_reset).cols()<<endl;
                    (*coarseMassMatrix_reset) = m_coarseP*(*coarseMassMatrix_reset)*m_coarseP.transpose();
                    
                    //        lambda can't capture member variable, so create a local one for lambda in ASSEMBLELIST
                    AssemblerEigenSparseMatrix<double> &coarseStiffnessMatrix_reset = m_coarseStiffnessMatrix_reset;
                    
                    //get stiffness matrix
                    ASSEMBLEMATINIT(coarseStiffnessMatrix_reset, coarse_world_reset.getNumQDotDOFs(), coarse_world_reset.getNumQDotDOFs());
                    ASSEMBLELIST(coarseStiffnessMatrix_reset, coarse_world_reset.getSystemList(), getStiffnessMatrix);
                    ASSEMBLELIST(coarseStiffnessMatrix_reset, coarse_world_reset.getForceList(), getStiffnessMatrix);
                    ASSEMBLEEND(coarseStiffnessMatrix_reset);
                    
                    
                    //constraint Projection
                    (*coarseStiffnessMatrix_reset) = m_coarseP*(*coarseStiffnessMatrix_reset)*m_coarseP.transpose();
                    
                    //Eigendecomposition for the coarse mesh reset
                    std::pair<Eigen::MatrixXx<double>, Eigen::VectorXx<double> > m_coarseUs_reset;
                    m_coarseUs_reset = generalizedEigenvalueProblemNegative(((*coarseStiffnessMatrix_reset)), (*m_coarseMassMatrix_reset), m_numModes, 0.00);
                    
                    // output for visual check
                    //                    unsigned int mode = 0;
                    //                    Eigen::VectorXd coarse_eig_def;
                    //                    for (mode = 0; mode < m_numModes; ++mode) {
                    //                        coarse_eig_def = (m_coarseP.transpose()*m_coarseUs_reset.first.col(mode)).transpose();
                    //                        idx = 0;
                    //                        // getGeometry().first is V
                    //                        Eigen::MatrixXd V_disp = std::get<0>(m_coarseWorld_reset.getSystemList().getStorage())[0]->getGeometry().first;
                    //                        for(unsigned int vertexId=0;  vertexId < std::get<0>(m_coarseWorld_reset.getSystemList().getStorage())[0]->getGeometry().first.rows(); ++vertexId) {
                    //
                    //                            V_disp(vertexId,0) += (1*coarse_eig_def(idx));
                    //                            idx++;
                    //                            V_disp(vertexId,1) += (1*coarse_eig_def(idx));
                    //                            idx++;
                    //                            V_disp(vertexId,2) += (1*coarse_eig_def(idx));
                    //                            idx++;
                    //                        }
                    
                    //                        Eigen::MatrixXi coarse_F = surftri(std::get<0>(m_coarseWorld_reset.getSystemList().getStorage())[0]->getGeometry().first, std::get<0>(m_coarseWorld_reset.getSystemList().getStorage())[0]->getGeometry().second);
                    //                        igl::writeOBJ("coarse_mesh_eigen_mode_reset" + std::to_string(mode) + ".obj",V_disp,coarse_F);
                    //
                    //                    }
                    //
                    for(int i = 0; i < m_numModes; ++i)
                    {
                        m_R(i) = m_fineUs_reset.second(i)/m_coarseUs_reset.second(i);
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
                        //#ifdef EDWIN_DEBUG
                        cout<<"calculated restart ratio:\n";
                        std::cout<<m_R(i)<<std::endl;
                        //#endif
                        
                    }
                    ratio_calculated = true;
                
                
            }
            else if (ratio_recalculation_switch == 3)
            {
                
                if(step_number%1 == 0) // update once every few frames
                {
                    Eigen::VectorXd eigv;
                    if(Eigen::loadMarketVector(eigv, "ExampleEigenFitRIMskater_5skater_3 100000 0.15 8 0 150 0 20 1 0.01 0 0 -0.01/fineeigenvalues" + std::to_string(step_number) + ".mtx"))
                    {
                        cout<<"loading eigenvalues from fine mesh at frame "<<step_number<<endl;
                        
                        
                        for(int i = 0; i < m_numModes; ++i)
                        {
                            m_R(i) = eigv(i)/m_coarseUs.second(i);
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
                            //#ifdef EDWIN_DEBUG
                            std::cout<<m_R(i)<<std::endl;
                            //#endif
                            
                        }
                        
                        ratio_calculated = true;
                    }
                    else{
                        cout<<"can't load eigenvalues from fine mesh"<<endl;
                        exit(1);
                    }
                }
                
                
            }
            else if (ratio_recalculation_switch == 4)
            {
                
                if(step_number%1 == 0) // update once every few frames
                {
                    
                    Eigen::MatrixXd fine_V;
                    Eigen::MatrixXd fine_F;
                    if(igl::readOBJ("FineCompSimExampleEigenFitRIMskater_3 100000 0.15 8 0 300 0 0 1 0.01 0 0 -0.01/pos" + std::to_string(step_number) + ".obj",fine_V, fine_F))
                    {
                        cout<<"loading fine mesh at frame "<<step_number<<endl;
                        
                        Eigen::VectorXd fine_displacement(3*fine_V.rows());
                        fine_displacement.setZero();
                        //                        fine_eig_def = (m_fineP.transpose()*m_fineUs_reset.first.col(mode)).transpose();
                        int idx = 0;
                        // getGeometry().first is V
                        Eigen::MatrixXd fine_V_disp0 = std::get<0>(m_fineWorld.getSystemList().getStorage())[0]->getGeometry().first;
                        for(unsigned int vertexId=0;  vertexId < std::get<0>(m_fineWorld.getSystemList().getStorage())[0]->getGeometry().first.rows(); ++vertexId) {
                            
                            //                            fine_displacement(idx) = fine_V(vertexId,0) - fine_V_disp0(vertexId,0);
                            fine_displacement(idx) = fine_V(vertexId,0);
                            //                            fine_V_disp(vertexId,0) += (1*fine_eig_def(idx));
                            idx++;
                            
                            //                            fine_displacement(idx) = fine_V(vertexId,1) - fine_V_disp0(vertexId,1);
                            fine_displacement(idx) = fine_V(vertexId,1);
                            //                            fine_V_disp(vertexId,1) += (1*fine_eig_def(idx));
                            idx++;
                            
                            //                            fine_displacement(idx) = fine_V(vertexId,2) - fine_V_disp0(vertexId,2);
                            fine_displacement(idx) = fine_V(vertexId,2);
                            //                            fine_V_disp(vertexId,2) += (1*fine_eig_def(idx));
                            idx++;
                        }
                        
                        World<double, std::tuple<PhysicalSystemImpl *>,
                        std::tuple<ForceSpringFEMParticle<double> *, ForceParticlesGravity<double> *>,
                        std::tuple<ConstraintFixedPoint<double> *> > &world = m_fineWorld;
                        
                        Eigen::Map<Eigen::VectorXd> fine_q = mapStateEigen<0>(m_fineWorld);
                        
                        fine_q = fine_displacement;
                        //            double pd_fine_pos[world.getNumQDOFs()]; // doesn't work for MSVS
                        Eigen::Map<Eigen::VectorXd> eigen_fine_pos0(fine_pos0,world.getNumQDOFs());
                        
                        //                        Eigen::VectorXx<double> posFull;
                        //                        posFull = this->getFinePositionFull(q);
                        //
                        fine_q = fine_displacement - eigen_fine_pos0;
                        
                        //        lambda can't capture member variable, so create a local one for lambda in ASSEMBLELIST
                        AssemblerEigenSparseMatrix<double> &fineStiffnessMatrix = m_fineStiffnessMatrix;
                        
                        //            std::cout<<
                        
                        //get stiffness matrix
                        ASSEMBLEMATINIT(fineStiffnessMatrix, world.getNumQDotDOFs(), world.getNumQDotDOFs());
                        ASSEMBLELIST(fineStiffnessMatrix, world.getSystemList(), getStiffnessMatrix);
                        ASSEMBLELIST(fineStiffnessMatrix, world.getForceList(), getStiffnessMatrix);
                        ASSEMBLEEND(fineStiffnessMatrix);
                        
                        
                        //constraint Projection
                        (*fineStiffnessMatrix) = m_fineP*(*fineStiffnessMatrix)*m_fineP.transpose();
                        
                        
                        //Eigendecomposition for the embedded fine mesh
                        std::pair<Eigen::MatrixXx<double>, Eigen::VectorXx<double> > m_Us;
                        m_Us = generalizedEigenvalueProblemNegative(((*fineStiffnessMatrix)), (*m_fineMassMatrix), m_numModes, 0.00);
                        
                        fineEigMassProj = m_Us;
                        fineEig = m_Us;
                        fineEigMassProj.first = (*m_fineMassMatrix)*fineEigMassProj.first;
                        
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
                            //#ifdef EDWIN_DEBUG
                            std::cout<<m_R(i)<<std::endl;
                            //#endif
                            
                        }
                        
                        ratio_calculated = true;
                    }
                    else{
                        cout<<"can't load eigenvalues from fine mesh"<<endl;
                        exit(1);
                    }
                }
                
                
            }
            
            //            std::cout<<m_coarseUs.second<<std::endl;
            
        }
        Y = (*coarseMassMatrix)*m_coarseUs.first*(m_R-m_I).asDiagonal();
        Z =  (m_coarseUs.second.asDiagonal()*m_coarseUs.first.transpose()*(*coarseMassMatrix));
        //        Eigen::saveMarket(Y, "Y.dat");
        //        Eigen::saveMarket(Z, "Z.dat");
        
        // no error
        return 0;
    }
    
    
    Eigen::SparseMatrix<double> fixedPointProjectionMatrixCoarse(Eigen::VectorXi &indices) {
        
        std::vector<Eigen::Triplet<double> > triplets;
        Eigen::SparseMatrix<double> P;
        Eigen::VectorXi sortedIndices = indices;
        std::sort(sortedIndices.data(), sortedIndices.data()+indices.rows());
        
        //build a projection matrix P which projects fixed points out of a physical syste
        int fIndex = 0;
        
        //total number of DOFS in system
        
        unsigned int n = 3*this->getImpl().getV().rows();
        unsigned int m = n - 3*indices.rows();
        
        P.resize(m,n);
        
        //number of unconstrained DOFs
        unsigned int rowIndex =0;
        for(unsigned int vIndex = 0; vIndex < this->getImpl().getV().rows(); vIndex++) {
            
            while((vIndex < this->getImpl().getV().rows()) && (fIndex < sortedIndices.rows()) &&(vIndex == sortedIndices[fIndex])) {
                fIndex++;
                vIndex++;
            }
            
            if(vIndex == this->getImpl().getV().rows())
                break;
            
            //add triplet into matrix
            triplets.push_back(Eigen::Triplet<double>(this->getQ().getGlobalId() +rowIndex, this->getQ().getGlobalId() + 3*vIndex,1));
            triplets.push_back(Eigen::Triplet<double>(this->getQ().getGlobalId() +rowIndex+1, this->getQ().getGlobalId() + 3*vIndex+1, 1));
            triplets.push_back(Eigen::Triplet<double>(this->getQ().getGlobalId() +rowIndex+2, this->getQ().getGlobalId() + 3*vIndex+2, 1));
            
            rowIndex+=3;
        }
        
        P.setFromTriplets(triplets.begin(), triplets.end());
        
        //build the matrix and  return
        return P;
    }
    
    void writeSimpleMesh(const std::string mesh_file_name,
                         const Eigen::MatrixXd & V,
                         const Eigen::MatrixXi & T)
    {
        // copied from tetwild
        std::fstream f(mesh_file_name, std::ios::out);
        f.precision(std::numeric_limits<double>::digits10 + 1);
        f << "MeshVersionFormatted 1" << std::endl;
        f << "Dimension 3" << std::endl;
        
        f << "Vertices" << " " <<V.rows()  << std::endl;
        for (int i = 0; i < V.rows(); i++)
            f << V(i,0) << " " << V(i,1) << " " << V(i,2) << " " << 0 << std::endl;
        f << "Tetrahedra" << std::endl;
        f << T.rows() << std::endl;
        for (int i = 0; i < T.rows(); i++) {
            for (int j = 0; j < 4; j++)
                f << T(i, j) + 1 << " ";
            f << 0 << std::endl;
        }
        
        f << "End";
        f.close();
    }
    
    Eigen::MatrixXi surftri(const Eigen::MatrixXd & V,
                            const Eigen::MatrixXi & T)
    {
        // translated from distmesh
        Eigen::MatrixXi faces(T.rows()*4,3);
        faces.col(0) << T.col(0),T.col(0),T.col(0),T.col(1);
        faces.col(1) << T.col(1),T.col(1),T.col(2),T.col(2);
        faces.col(2) << T.col(2),T.col(3),T.col(3),T.col(3);
        
        // the fourth vertex of the tet that's not on the surface
        Eigen::VectorXi node4(T.rows()*4);
        node4 << T.col(3), T.col(2), T.col(1), T.col(0);
        
        Eigen::MatrixXi facesSorted(T.rows()*4,3);
        Eigen::VectorXi SortedInd(T.rows()*4);
        igl::sortrows(faces,true,facesSorted,SortedInd);
        Eigen::MatrixXi C;
        Eigen::VectorXi IA;
        Eigen::VectorXi IC;
        igl::unique_rows(faces,C,IA,IC);
        
        Eigen::VectorXi count(IC.maxCoeff());
        Eigen::VectorXi histbin(IC.maxCoeff());
        //        std::cout<<IC.maxCoeff()<<std::endl;
        
        for (int ind = 0; ind < IC.maxCoeff(); ++ind) {
            histbin(ind) = ind;
        }
        //        std::cout<<histbin<<std::endl;
        //        std::cout<<IC<<std::endl;
        Eigen::VectorXi foo(IC.size());
        igl::histc(IC,histbin,count,foo);
        //        std::cout<<count<<std::endl;
        std::vector<int> nonDuplicatedFacesInd;
        int oneCount = 0;
        for (int ind = 0; ind < count.size(); ++ind) {
            if (count(ind)==1) {
                ++oneCount;
                nonDuplicatedFacesInd.push_back(ind);
            }
        }
        
        Eigen::MatrixXi nonDuplicatedFaces(nonDuplicatedFacesInd.size(),3);
        Eigen::VectorXi nonDuplicatedNode4(nonDuplicatedFacesInd.size());
        for(std::vector<int>::iterator it = nonDuplicatedFacesInd.begin(); it != nonDuplicatedFacesInd.end(); ++it)
        {
            nonDuplicatedFaces.row(it-nonDuplicatedFacesInd.begin()) << faces.row(IA(*it));
            nonDuplicatedNode4(it-nonDuplicatedFacesInd.begin()) = node4(IA(*it));
        }
        
        // use the three vectors of the tet to determine the orientation
        Eigen::Vector3d v1;
        Eigen::Vector3d v2;
        Eigen::Vector3d v3;
        int tempV;
        for (int ind = 0; ind < nonDuplicatedFaces.rows(); ++ind) {
            v1 = V.row(nonDuplicatedFaces(ind,1)) - V.row(nonDuplicatedFaces(ind,0));
            v2 = V.row(nonDuplicatedFaces(ind,2)) - V.row(nonDuplicatedFaces(ind,0));
            v3 = V.row(nonDuplicatedNode4(ind)) - V.row(nonDuplicatedFaces(ind,0));
            if ((v1.cross(v2)).dot(v3) > 0) {
                
                tempV = nonDuplicatedFaces(ind,1);
                nonDuplicatedFaces(ind,1) = nonDuplicatedFaces(ind,2);
                nonDuplicatedFaces(ind,2) = tempV;
            }
        }
        
        return nonDuplicatedFaces;
    }
    
    //per vertex accessors. takes the state of the coarse mesh
    inline Eigen::Vector3x<double> getFinePosition(const State<double> &state, unsigned int vertexId) const {
        return m_Vf.row(vertexId).transpose() + m_N.block(3*vertexId, 0, 3, m_N.cols())*(*this).getImpl().getElement(m_elements[vertexId])->q(state);
    }
    
    
    inline Eigen::VectorXx<double> getFinePositionFull(const Eigen::VectorXd V) const {
        //            std::cout<<m_Vf.rows()<<std::endl;
        //            std::cout<<m_Vf.cols()<<std::endl;
        //            std::cout<<V.rows()<<std::endl;
        //            std::cout<<V.cols()<<std::endl;
        //            std::cout<<((*N)*V).rows()<<std::endl;
        //            std::cout<<((*N)*V).cols()<<std::endl;
        Eigen::Map<Eigen::VectorXd> eigen_fine_pos0(fine_pos0,m_fineWorld.getNumQDOFs());
        
        return eigen_fine_pos0 + (*N)*V;
    }
    
    
    inline World<double, std::tuple<PhysicalSystemImpl *>,
    std::tuple<ForceSpringFEMParticle<double> *, ForceParticlesGravity<double> *>,
    std::tuple<ConstraintFixedPoint<double> *> > & getFineWorld(){ return m_fineWorld;}
    
    
    std::pair<Eigen::MatrixXx<double>, Eigen::VectorXx<double> > coarseEigMassProj;
    std::pair<Eigen::MatrixXx<double>, Eigen::VectorXx<double> > fineEigMassProj;
    std::pair<Eigen::MatrixXx<double>, Eigen::VectorXx<double> > coarseEig;
    std::pair<Eigen::MatrixXx<double>, Eigen::VectorXx<double> > fineEig;
    
    Eigen::MatrixXd fforthogonal;
    Eigen::MatrixXd cforthogonal;
    Eigen::MatrixXd fforthogonal2;
    Eigen::MatrixXd cforthogonal2;
    Eigen::MatrixXd fforthogonal3;
    Eigen::MatrixXd cforthogonal3;
    Eigen::MatrixXd fforthogonal4;
    Eigen::MatrixXd cforthogonal4;
    Eigen::MatrixXd fforthogonal5;
    Eigen::MatrixXd cforthogonal5;
    Eigen::MatrixXd fforthogonal6;
    Eigen::MatrixXd cforthogonal6;
    
    Eigen::SparseMatrix<double> m_fineP;
    Eigen::SparseMatrix<double> m_coarseP;
    
    //        Eigen::MatrixXd coarse_V_disp_p;
    
    AssemblerEigenSparseMatrix<double> m_coarseMassMatrix;
    AssemblerEigenSparseMatrix<double> m_fineMassMatrix;
    
    AssemblerEigenVector<double> m_fineforceVector;
    AssemblerEigenVector<double> m_finefExt;
    
    AssemblerEigenVector<double> m_forceVector;
    
    AssemblerEigenSparseMatrix<double> N;
    
    SolverPardiso<Eigen::SparseMatrix<double, Eigen::RowMajor> > m_pardiso_test;
    Eigen::VectorXd minvf;
    Eigen::VectorXd minvfCP;
    
    // rest state of fine q
    double* fine_pos0  = NULL;
    // rest state of coarse q
    double* coarse_pos0 = NULL;
    
    double* fine_q_transfered = NULL;
    
    bool haus = false;
    
    int flag = 0;
    
    // TODO: parameters for rayleigh damping. should not be here...
    double a;
    double b;
    
    int step_number;
    bool ratio_calculated;
protected:
    
    World<double, std::tuple<PhysicalSystemImpl *>,
    std::tuple<ForceSpringFEMParticle<double> *, ForceParticlesGravity<double> *>,
    std::tuple<ConstraintFixedPoint<double> *> > m_fineWorld;
    
    World<double, std::tuple<PhysicalSystemImpl *>,
    std::tuple<ForceSpringFEMParticle<double> *, ForceParticlesGravity<double> *>,
    std::tuple<ConstraintFixedPoint<double> *> > m_coarseWorld;
    
    
    AssemblerEigenSparseMatrix<double> m_fineStiffnessMatrix;
    
    AssemblerEigenSparseMatrix<double> m_coarseStiffnessMatrix;
    
    Eigen::VectorXi m_fineMovingVerts;
    Eigen::VectorXi m_fineFixedVerts;
    
    Eigen::VectorXi m_coarseMovingVerts;
    Eigen::VectorXi m_coarseFixedVerts;
    
    
    double youngs;
    double poisson;
    
    int m_constraintDir;
    double m_constraintTol;
    
    //num modes to correct
    unsigned int m_numModes;
    //    unsigned int m_numToCorrect;
    
    //Ratios diagonal matrix, stored as vector
    Eigen::VectorXd m_R;
    Eigen::VectorXd m_I;
    
    
    
    Eigen::MatrixXx<double> m_Vf;
    Eigen::MatrixXi m_Ff;
    Eigen::MatrixXd m_N;
    //m_elements[i] is a n-vector that stores the index of the element containing the ith vertex in the embedded mesh
    Eigen::VectorXi m_elements;
    
    State<double> restFineState;
    
    
    bool ratio_recalculation_flag;
    int ratio_recalculation_switch;
    
    
    unsigned int constraint_switch;
    unsigned int m_numConstraints;
private:
    
};

#endif /* EigenFit_h */

