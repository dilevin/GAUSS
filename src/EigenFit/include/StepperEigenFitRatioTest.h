//
//  StepperEigenFitRatioTest.h
//  Gauss
//
//  Created by Edwin Chen on 2018-05-15.
//
//

#ifndef StepperEigenFitRatioTest_h
#define StepperEigenFitRatioTest_h



#include <World.h>
#include <Assembler.h>
#include <TimeStepper.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <UtilitiesEigen.h>
#include <UtilitiesMATLAB.h>
#include <Eigen/SparseCholesky>
#include <SolverPardiso.h>
#include <EigenFit.h>
#include <limits>
#include <igl/boundary_facets.h>


//TODO Solver Interface
namespace Gauss {
    
    template<typename DataType, typename MatrixAssembler, typename VectorAssembler>
    class StepperEigenFitRatioTestImpl
    {
    public:
        
        template<typename Matrix>
        StepperEigenFitRatioTestImpl(Matrix &P, unsigned int numModes, std::string cmeshname, std::string fmeshname) {
            
            m_numModes = (numModes);
            
            //            std::cout<<m_P.rows()<<std::endl;
            m_P = P;
            m_factored = false;
            // refactor for every solve
            m_refactor = true;
            
            // init residual
            res = std::numeric_limits<double>::infinity();
            
            it_outer = 0;
            it_inner = 0;
            // constants from Nocedal and Wright
            step_size = 1;
            c1 = 1e-4;
            c2 = 0.9;
            
            frame = 1;
            
            m_fmeshname = fmeshname;
            m_cmeshname = cmeshname;
        }
        
        StepperEigenFitRatioTestImpl(const StepperEigenFitRatioTestImpl &toCopy) {
            
        }
        
        ~StepperEigenFitRatioTestImpl() { }
        
        //Methods
        //init() //initial conditions will be set at the begining
        template<typename World>
        void step(World &world, double dt, double t);
        
        inline typename VectorAssembler::MatrixType & getLagrangeMultipliers() { return m_lagrangeMultipliers; }
        
        double* rhs  = NULL;
        double* v_old  = NULL;
        double* v_temp = NULL;
        double* q_old = NULL;
        
    protected:
        
        //num modes to correct
        unsigned int m_numModes;
        
        //        //Ratios diagonal matrix, stored as vector
        Eigen::VectorXd m_R;
        
        
        //Subspace Eigenvectors and eigenvalues from the embedded fine mesh
        std::pair<Eigen::MatrixXx<DataType>, Eigen::VectorXx<DataType> > m_fineUs;
        //Subspace Eigenvectors and eigenvalues from this coarse mesh
        std::pair<Eigen::MatrixXx<DataType>, Eigen::VectorXx<DataType> > m_coarseUs;
        // for SMW
        Eigen::MatrixXd Y;
        Eigen::MatrixXd Z;
        
        MatrixAssembler m_massMatrix;
        MatrixAssembler m_stiffnessMatrix;
        MatrixAssembler m_fineMassMatrix;
        MatrixAssembler m_fineStiffnessMatrix;
        VectorAssembler m_forceVector;
        VectorAssembler m_fExt;
        
        //        // for calculating the residual. ugly
        //        MatrixAssembler m_massMatrixNew;
        //        MatrixAssembler m_stiffnessMatrixNew;
        //        VectorAssembler m_forceVectorNew;
        //        VectorAssembler m_fExtNew;
        
        
        Eigen::SparseMatrix<DataType> m_P;
        
        //storage for lagrange multipliers
        typename VectorAssembler::MatrixType m_lagrangeMultipliers;
        
        bool m_factored, m_refactor;
        
        // iteration counter
        int it_outer, it_inner;
        // residual
        double res, res_old, step_size, c1, c2;
        
        int frame;
        
        std::string m_fmeshname, m_cmeshname;
        //        Eigen::VectorXd res;
        
#ifdef GAUSS_PARDISO
        
        SolverPardiso<Eigen::SparseMatrix<DataType, Eigen::RowMajor> > m_pardiso;
        SolverPardiso<Eigen::SparseMatrix<DataType, Eigen::RowMajor> > m_pardiso_test;
        SolverPardiso<Eigen::SparseMatrix<DataType, Eigen::RowMajor> > m_pardiso_mass;
        //        SolverPardiso<Eigen::SparseMatrix<DataType, Eigen::RowMajor> > m_pardiso_res;
#else
#endif
        
    private:
    };
}

template<typename DataType, typename MatrixAssembler, typename VectorAssembler>
template<typename World>
void StepperEigenFitRatioTestImpl<DataType, MatrixAssembler, VectorAssembler>::step(World &world, double dt, double t) {
    
    //First two lines work around the fact that C++11 lambda can't directly capture a member variable.
    VectorAssembler &forceVector = m_forceVector;
    VectorAssembler &fExt = m_fExt;
    
    // init rhs
    rhs = new double[m_P.rows()];
    Eigen::Map<Eigen::VectorXd> eigen_rhs(rhs,m_P.rows());
    
    // velocity from previous step (for calculating residual)
    v_old = new double[world.getNumQDotDOFs()];
    Eigen::Map<Eigen::VectorXd> eigen_v_old(v_old,world.getNumQDotDOFs());
    
    // velocity from previous step (for calculating residual)
    v_temp = new double[world.getNumQDotDOFs()];
    Eigen::Map<Eigen::VectorXd> eigen_v_temp(v_temp,world.getNumQDotDOFs());
    
    // velocity from previous step (for calculating residual)
    q_old = new double[world.getNumQDotDOFs()];
    Eigen::Map<Eigen::VectorXd> eigen_q_old(q_old,world.getNumQDotDOFs());
    
    //Grab the state
    //        Eigen::VectorXd qDotOld = qDot; // deep copy
    
    
    Eigen::Map<Eigen::VectorXd> fine_q =   mapStateEigen<0>(static_cast<EigenFit*>(std::get<0>(world.getSystemList().getStorage())[0])->getFineWorld());
    
    
    std::string qfileName(m_fmeshname + "_def" + std::to_string(frame) + ".mtx");
    Eigen::VectorXd  tempv;
    Eigen::loadMarketVector(tempv,qfileName);
    
//    std::cout<< tempv.rows() << std::endl;
    
//    std::cout<<"original state size "<<q.rows()<<"\nloaded state size "<<tempv.rows()<<"\n";
//    std::cout<<"fine state size "<<fine_q.rows()<<"\nloaded state size "<<tempv.rows()<<"\n";
    fine_q = tempv.head(fine_q.rows());
    
    static_cast<EigenFit*>(std::get<0>(world.getSystemList().getStorage())[0])->calculateFineMesh();
    
//    std::cout<<(static_cast<EigenFit*>(std::get<0>(world.getSystemList().getStorage())[0])->fineEig).second<<std::endl;
//      igl::writeOBJ("finesurfpos.obj",Vf_disp,surfFf);
//
//    std::cout<<"N size "<<(*(static_cast<EigenFit*>(std::get<0>(world.getSystemList().getStorage())[0])->N)).rows()<<" "<<(*(static_cast<EigenFit*>(std::get<0>(world.getSystemList().getStorage())[0])->N)).cols()<<"\n";
    
    Eigen::Map<Eigen::VectorXd> q = mapStateEigen<0>(world);
    
//
    std::string cqfileName(m_cmeshname + "_def" + std::to_string(frame) + ".mtx");
    Eigen::VectorXd  ctempv;
    Eigen::loadMarketVector(ctempv,cqfileName);
    q = ctempv.head(q.rows());
    
//    q = (*(static_cast<EigenFit*>(std::get<0>(world.getSystemList().getStorage())[0])->N)).transpose() * fine_q;
//


    //First two lines work around the fact that C++11 lambda can't directly capture a member variable.
    MatrixAssembler &massMatrix = m_massMatrix;
    MatrixAssembler &stiffnessMatrix = m_stiffnessMatrix;

    //get mass matrix
    ASSEMBLEMATINIT(massMatrix, world.getNumQDotDOFs(), world.getNumQDotDOFs());
    ASSEMBLELIST(massMatrix, world.getSystemList(), getMassMatrix);
    ASSEMBLEEND(massMatrix);


    //get stiffness matrix
    ASSEMBLEMATINIT(stiffnessMatrix, world.getNumQDotDOFs(), world.getNumQDotDOFs());
    ASSEMBLELIST(stiffnessMatrix, world.getSystemList(), getStiffnessMatrix);
    ASSEMBLELIST(stiffnessMatrix, world.getForceList(), getStiffnessMatrix);
    ASSEMBLEEND(stiffnessMatrix);

    //constraint Projection
    (*massMatrix) = m_P*(*massMatrix)*m_P.transpose();
    (*stiffnessMatrix) = m_P*(*stiffnessMatrix)*m_P.transpose();

    std::pair<Eigen::MatrixXx<double>, Eigen::VectorXx<double> > m_coarseUs;
    // matrices passed in already eliminated the constraints
    // Eigendecomposition for the coarse mesh
    m_coarseUs = generalizedEigenvalueProblem((*stiffnessMatrix), (*massMatrix), m_numModes,0.000);
    
    
    

    for (int i = 0; i < m_numModes; i ++) {
        std::cout<<(static_cast<EigenFit*>(std::get<0>(world.getSystemList().getStorage())[0])->fineEig).second(i) / (m_coarseUs).second(i) << std::endl;
        
        
        
    }
    std::cout<<"-\n";
    
    Eigen::saveMarket((static_cast<EigenFit*>(std::get<0>(world.getSystemList().getStorage())[0])->fineEig).second, "fe" + std::to_string(frame) + ".dat");
    
    Eigen::saveMarket( (m_coarseUs).second, "ce" + std::to_string(frame) + ".dat");
    
    // embedded fine pos
    fine_q = (*(static_cast<EigenFit*>(std::get<0>(world.getSystemList().getStorage())[0])->N)) * q;
    static_cast<EigenFit*>(std::get<0>(world.getSystemList().getStorage())[0])->calculateFineMesh();
    
    for (int i = 0; i < m_numModes; i ++) {
        std::cout<<(static_cast<EigenFit*>(std::get<0>(world.getSystemList().getStorage())[0])->fineEig).second(i) / (m_coarseUs).second(i) << std::endl;
        
        
        
    }
    std::cout<<"--\n";
    
    Eigen::saveMarket((static_cast<EigenFit*>(std::get<0>(world.getSystemList().getStorage())[0])->fineEig).second, "fee" + std::to_string(frame) + ".dat");
    
//
//    for (int i = 0; i < m_numModes; i ++) {
//        std::cout<<(m_coarseUs).second(i) << std::endl;
//    }
//    std::cout<<"--\n";
    
    
    frame++;
    
}

template<typename DataType, typename MatrixAssembler, typename VectorAssembler>
using StepperEigenFitRatioTest = TimeStepper<DataType, StepperEigenFitRatioTestImpl<DataType, MatrixAssembler, VectorAssembler> >;




#endif /* StepperEigenFitRatioTest_h */
