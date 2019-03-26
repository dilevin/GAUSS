//
//  UtilitiesEigen.h
//  Gauss
//
//  Created by David Levin on 2/11/17.
//
//

#ifndef UtilitiesEigen_h
#define UtilitiesEigen_h

#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/src/SparseCore/SparseSelfAdjointView.h>
#include <Utilities.h>
#include <World.h>

//Modal Analysis using Spectra
#include <GenEigsComplexShiftSolver.h>
#include <SymGEigsSolver.h>
#include <GenEigsRealShiftSolver.h>
#include <MatOp/SparseGenMatProd.h>
#include <MatOp/SparseCholesky.h>
#include <MatOp/SparseSymShiftSolve.h>
#include <SymGEigsSolver.h>
#include <stdexcept>

//some useful types
namespace Eigen {
    template<typename DataType>
    using Vector3x = Eigen::Matrix<DataType, 3,1>;
    
    template<typename DataType>
    using Vector6x = Eigen::Matrix<DataType, 6,1>;
    
    template<typename DataType>
    using VectorXx = Eigen::Matrix<DataType, Eigen::Dynamic, 1>;
    
    template<typename DataType>
    using MatrixXx = Eigen::Matrix<DataType, Eigen::Dynamic, Eigen::Dynamic>;
    
    template<typename DataType>
    using Matrix33x = Eigen::Matrix<DataType, 3, 3>;
    
    //Comes up a lot in constitutive models
    template<typename DataType>
    using Matrix66x = Eigen::Matrix<DataType, 6,6>;
    
    template<typename DataType>
    using Tensor333x = std::array<std::array<Eigen::Vector3x<DataType>,3>, 3>;
    
    template<typename DataType>
    using Tensor3333x = std::array<std::array<Eigen::Matrix33x<DataType>,3>, 3>;
    
    //useful maps
    template<typename DataType>
    using Map3x = Eigen::Map<Vector3x<DataType> >;
    
    
    //Derivative of 3x3 SVD
    //takes in the svd of a 3x3 matrix (U,S,V s.t A = U*S*V^T) where U,V are orthornomal, S is diagonal stored as a vector.
    //The derivative uses fast SVD code by Eftychois Sifakis via libigl which only works corectly for single precision floating point values
    //Inputs
    //
    //Output
    template<typename DataType,
             typename DerivedU,
             typename DerivedV,
             typename DerivedS>
    inline void dSVD(Tensor3333x<DataType> &dU,
                     Tensor333x<DataType>  &dS,
                     Tensor3333x<DataType> &dV,
                     const Eigen::MatrixBase<DerivedU> &U,
                     const Eigen::MatrixBase<DerivedS> &S,
                     const Eigen::MatrixBase<DerivedV> &V) {
        
        Eigen::Matrix33x<typename DerivedU::Scalar> UVT, tmp;
        
        Eigen::Matrix33x<typename DerivedU::Scalar> lambda;
    
        //crappy hack for now
        typename DerivedU::Scalar tol = 1e-6;
        typename DerivedU::Scalar w01, w02, w12;
        typename DerivedU::Scalar d01, d02, d12;
        
        d01 = S(1)*S(1)-S(0)*S(0);
        d02 = S(2)*S(2)-S(0)*S(0);
        d12 = S(2)*S(2)-S(1)*S(1);
        
        //corresponds to conservative solution --- if singularity is detected no angular velocity
        d01 = static_cast<typename DerivedU::Scalar>(1.0)/(std::abs(d01) < tol ? std::numeric_limits<double>::infinity() : d01);
        d02 = static_cast<typename DerivedU::Scalar>(1.0)/(std::abs(d02) < tol ? std::numeric_limits<double>::infinity() : d02);
        d12 = static_cast<typename DerivedU::Scalar>(1.0)/(std::abs(d12) < tol ? std::numeric_limits<double>::infinity() : d12);
        
        for(unsigned int r=0; r<3; ++r) {
            for(unsigned int s =0; s <3; ++s) {
                
                UVT = U.row(r).transpose()*V.row(s);
                
                //Compute dS
                dS[r][s] = UVT.diagonal();
                
                UVT -= dS[r][s].asDiagonal();
                
                tmp  = S.asDiagonal()*UVT + UVT.transpose()*S.asDiagonal();
                w01 = tmp(0,1)*d01;
                w02 = tmp(0,2)*d02;
                w12 = tmp(1,2)*d12;
                tmp << 0, w01, w02,
                       -w01, 0, w12,
                       -w02, -w12, 0;
                
                dV[r][s] = V*tmp;
                
                tmp = UVT*S.asDiagonal() + S.asDiagonal()*UVT.transpose();
                w01 = tmp(0,1)*d01;
                w02 = tmp(0,2)*d02;
                w12 = tmp(1,2)*d12;
                tmp << 0, w01, w02,
                -w01, 0, w12,
                -w02, -w12, 0;
                
                dU[r][s] = U*tmp;
                
                //compute dV
                //dU[r][s] = lambda.fullPivLu().solve(S.asDiagonal()*UVT - UVT.transpose()*S.asDiagonal());
                //dU[r][s] = V.transpose()*dV[r][s];
                
            }
        }
    }
    
}

namespace Gauss {
    //state ptr direct to eigen map for a single property (position or velocity)
    template<unsigned int Property, typename World>
    Eigen::Map<Eigen::VectorXx<typename World::Scalar> > mapStateEigen(World &world) {
        std::tuple<typename World::Scalar *, unsigned int> ptr = world.getState().template getStatePtr<Property>();
        return Eigen::Map<Eigen::VectorXx<typename World::Scalar> >(std::get<0>(ptr), std::get<1>(ptr));
    }
    
    template<unsigned int Property, typename DataType>
    Eigen::Map<Eigen::VectorXx<typename Gauss::State<DataType>::StateDataType> > mapStateEigen(Gauss::State<DataType> &state) {
        std::tuple<typename Gauss::State<DataType>::StateDataType *, unsigned int> ptr = state.template getStatePtr<Property>();
        return Eigen::Map<Eigen::VectorXx<typename Gauss::State<DataType>::StateDataType> >(std::get<0>(ptr), std::get<1>(ptr));
    }
    
    template<unsigned int Property, typename DataType>
    Eigen::Map<Eigen::VectorXx<typename Gauss::State<DataType>::StateDataType> > mapStateEigen(const Gauss::State<DataType> &state) {
        std::tuple<typename Gauss::State<DataType>::StateDataType *, unsigned int> ptr = state.template getStatePtr<Property>();
        return Eigen::Map<Eigen::VectorXx<typename Gauss::State<DataType>::StateDataType> >(std::get<0>(ptr), std::get<1>(ptr));
    }
    
    //state ptr for the whole thing
    template<typename World>
    Eigen::Map<Eigen::VectorXx<typename World::Scalar> > mapStateEigen(World &world) {
        std::tuple<typename World::Scalar *, unsigned int> ptr = world.getState().getStatePtr();
        return Eigen::Map<Eigen::VectorXx<typename World::Scalar> >(std::get<0>(ptr), std::get<1>(ptr));
    }
    
    template<typename DOF, typename DataType, typename ...SystemTypes, typename ...ForceTypes, typename ...ConstraintTypes>
    inline Eigen::Map<Eigen::VectorXd> mapDOFEigen(DOF &dof, World<DataType,
                                                   std::tuple<SystemTypes...>,
                                                   std::tuple<ForceTypes...>,
                                                   std::tuple<ConstraintTypes...> > &world) {
        std::tuple<double *, unsigned int> qPtr = dof.getPtr(world.getState());
        //set position DOF and check
        return Eigen::Map<Eigen::VectorXd>(std::get<0>(qPtr), dof.getNumScalarDOF());
        
    }
    
    template<typename DOF, typename DataType>
    inline Eigen::Map<Eigen::VectorXx<DataType> > mapDOFEigen(DOF &dof, const State<DataType> &state) {
        std::tuple<DataType *, unsigned int> qPtr = dof.getPtr(state);
        //set position DOF and check
        return Eigen::Map<Eigen::VectorXx<DataType> >(std::get<0>(qPtr), dof.getNumScalarDOF());
    }
    
    //functor for getting position of a DOF
    template <typename DataType, typename DOF>
    class PositionEigen {
    public:
        inline PositionEigen(DOF *dof=nullptr, Eigen::Vector3x<DataType> p = Eigen::Vector3x<DataType>()) { m_dof = dof; m_p = p; }
        inline Eigen::Vector3x<DataType> operator()(const State<DataType> &state) const { return m_p+ mapDOFEigen(*m_dof, state); }
        inline DOF * getDOF() { return m_dof; }
    protected:
        DOF *m_dof;
        Eigen::Vector3x<DataType> m_p;
    };
}
//Define a new spectra shift and invert for "mass" shifting the generalized eigenproblem (cut and paste from Spectra's SymShiftSolve
namespace Spectra {
    
    
    ///
    /// \ingroup MatOp
    ///
    /// This class defines the shift-solve operation on a sparse real symmetric matrix \f$A\f$,
    /// i.e., calculating \f$y=(A-\sigma I)^{-1}x\f$ for any real \f$\sigma\f$ and
    /// vector \f$x\f$. It is mainly used in the SymEigsShiftSolver eigen solver.
    ///
    template <typename Scalar, int Uplo = Eigen::Lower, int Flags = 0, typename StorageIndex = int>
    class SparseSymMassShiftSolve
    {
    private:
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
        typedef Eigen::Map<const Vector> MapConstVec;
        typedef Eigen::Map<Vector> MapVec;
        typedef Eigen::SparseMatrix<Scalar, Flags, StorageIndex> SparseMatrix;
        
        const SparseMatrix m_stiffnessMat, m_massMat;
        const int m_n;
        Eigen::SimplicialLDLT<SparseMatrix, Uplo> m_solver;
        
        
    public:
        ///
        /// Constructor to create the matrix operation object.
        ///
        /// \param mat_ An **Eigen** sparse matrix object, whose type is
        /// `Eigen::SparseMatrix<Scalar, ...>`.
        ///
        SparseSymMassShiftSolve(const SparseMatrix& stiffnessMat_, const SparseMatrix &massMat_) :
        m_stiffnessMat(stiffnessMat_), m_massMat(massMat_),
        m_n(massMat_.rows())
        {
            if(stiffnessMat_.rows() != stiffnessMat_.cols() || massMat_.rows() != massMat_.cols())
                throw std::invalid_argument("SparseSymMassShiftSolve: matrices must be square");
            
            if(stiffnessMat_.rows() != massMat_.rows()) {
                throw std::invalid_argument("SparseSymMassShiftSolve: matrices must be the same size");
            }
        }
        
        ///
        /// Return the number of rows of the underlying matrix.
        ///
        int rows() const { return m_n; }
        ///
        /// Return the number of columns of the underlying matrix.
        ///
        int cols() const { return m_n; }
        
        ///
        /// Set the real shift \f$\sigma\f$.
        ///
        void set_shift(Scalar sigma)
        {
            //sadly there's no simple setShit here
            m_solver.compute(m_stiffnessMat+sigma*m_massMat);
            
            if(m_solver.info()!=Eigen::Success) {
                std::cout<<"Mass Shift: decomposition failed \n";
                exit(1);
            }
        }
        
        ///
        /// Perform the shift-solve operation \f$y=(A-\sigma I)^{-1}x\f$.
        ///
        /// \param x_in  Pointer to the \f$x\f$ vector.
        /// \param y_out Pointer to the \f$y\f$ vector.
        ///
        // y_out = inv(A - sigma * I) * x_in
        void perform_op(const Scalar* x_in, Scalar* y_out) const
        {
            MapConstVec x(x_in,  m_n);
            MapVec      y(y_out, m_n);
            y.noalias() = m_solver.solve(m_massMat*x);
        }
    };
    
    
} // namespace Spectra

//solve sparse generalized eigenvalue problem using spectra
//solve the gevp Ax = lambda*Bx
template<typename DataType, int Flags, typename Indices>
auto generalizedEigenvalueProblem(const Eigen::SparseMatrix<DataType, Flags, Indices> &A,
                                  const Eigen::SparseMatrix<DataType, Flags,Indices> &B,
                                  unsigned int numVecs) {
    
    //Spectra seems to freak out if you use row storage, this copy just ensures everything is setup the way the solver likes
    Eigen::SparseMatrix<DataType> K = A;
    Eigen::SparseMatrix<DataType> M = B;
    
    
    Spectra::SparseSymMatProd<DataType> Aop(K);
    Spectra::SparseCholesky<DataType>   Bop(M);
    
    //Spectra::SparseSymShiftSolve<DataType> Aop(K);
    
    //Aop.set_shift(1e-3);
    
    // Construct eigen solver object, requesting the smallest three eigenvalues
    Spectra::SymGEigsSolver<DataType, Spectra::SMALLEST_MAGN, Spectra::SparseSymMatProd<DataType>, Spectra::SparseCholesky<DataType>, Spectra::GEIGS_CHOLESKY > eigs(&Aop, &Bop, numVecs, 5*numVecs);
    
    
    // Initialize and compute
    eigs.init();
    //int nconv = eigs.compute();
    
    // Retrieve results
    if(eigs.info() == Spectra::SUCCESSFUL) {
        
        return std::make_pair(eigs.eigenvectors(), eigs.eigenvalues());
    } else {
        std::cout<<"Failure: "<<eigs.info()<<"\n";
        exit(1);
        return std::make_pair(eigs.eigenvectors(), eigs.eigenvalues());
    }
    
}

//use shift and invert to find Eigenvalues near the shift
template<typename DataType, int Flags, typename Indices>
auto generalizedEigenvalueProblem(const Eigen::SparseMatrix<DataType, Flags, Indices> &A,
                                  const Eigen::SparseMatrix<DataType, Flags,Indices> &B,
                                  unsigned int numVecs, DataType shift) {
    
    //Spectra seems to freak out if you use row storage, this copy just ensures everything is setup the way the solver likes
    Eigen::SparseMatrix<DataType> K = -A + shift*B;
    Eigen::SparseMatrix<DataType> M = B;
    
    //Spectra::SparseSymMassShiftSolve<DataType> Aop(K, M);
    //Aop.set_shift(shift);
    Spectra::SparseSymMatProd<DataType> Aop(M);
    Spectra::SparseCholesky<DataType> Bop(K);
    
    // Construct eigen solver object, requesting the smallest three eigenvalues
    Spectra::SymGEigsSolver<DataType, Spectra::LARGEST_MAGN, Spectra::SparseSymMatProd<DataType>, Spectra::SparseCholesky<DataType>, Spectra::GEIGS_CHOLESKY > eigs(&Aop, &Bop, numVecs, 5*numVecs);
    
    
    //Spectra::GenEigsShiftSolver<double, Spectra::LARGEST_MAGN,  Spectra::SparseSymMassShiftSolve<DataType> > eigs(&Aop, numVecs, 5*numVecs, shift);
    
    // Initialize and compute
    eigs.init();
    eigs.compute();
    Eigen::VectorXx<DataType> eigsCorrected;
    Eigen::MatrixXx<DataType> evsCorrected; //magnitude of eigenvectors can be wrong in this formulation
    eigsCorrected.resize(eigs.eigenvalues().rows());
    evsCorrected.resize(eigs.eigenvectors().rows(), eigs.eigenvectors().cols());

    
    // Retrieve results
    evsCorrected = eigs.eigenvectors();
    
    if(eigs.info() == Spectra::SUCCESSFUL) {
        //correct eigenvalues
        for(unsigned int ii=0; ii<eigs.eigenvalues().rows(); ++ii) {
            eigsCorrected[ii] = -(static_cast<DataType>(1)/(eigs.eigenvalues()[ii]) + shift);
            evsCorrected.col(ii) /= sqrt(eigs.eigenvectors().col(ii).transpose()*M*eigs.eigenvectors().col(ii));
        }
        
        return std::make_pair(evsCorrected, eigsCorrected);
    } else {
        std::cout<<"Failure: "<<eigs.info()<<"\n";
        exit(1);
        return std::make_pair(eigs.eigenvectors(), eigs.eigenvalues());
    }
    
}

#endif /* UtilitiesEigen_h */
