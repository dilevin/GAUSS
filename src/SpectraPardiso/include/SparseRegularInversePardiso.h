// Copyright (C) 2017-2018 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SPARSE_REGULAR_INVERSE_PARDISO_H
#define SPARSE_REGULAR_INVERSE_PARDISO_H

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
#include <stdexcept>
#include <UtilitiesEigen.h>
#include <SolverPardiso.h>

namespace Spectra {
    
    
    ///
    /// \ingroup MatOp
    ///
    /// This class defines matrix operations required by the generalized eigen solver
    /// in the regular inverse mode. For a sparse and positive definite matrix \f$B\f$,
    /// it implements the matrix-vector product \f$y=Bx\f$ and the linear equation
    /// solving operation \f$y=B^{-1}x\f$.
    ///
    /// This class is intended to be used with the SymGEigsSolver generalized eigen solver
    /// in the regular inverse mode.
    ///
    template <typename Scalar, int Uplo = Eigen::Lower, int Flags = 1, typename StorageIndex = int>
    class SparseRegularInversePardiso
    {
    private:
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
        typedef Eigen::Map<const Vector> MapConstVec;
        typedef Eigen::Map<Vector> MapVec;
        typedef Eigen::SparseMatrix<Scalar, Flags, StorageIndex> SparseMatrix;
        
        const int m_n;
        SparseMatrix& m_mat;
//        Eigen::ConjugateGradient<SparseMatrix> m_cg;
        
        mutable SolverPardiso<SparseMatrix> m_pardiso;


        
    public:
        
        ///
        /// Constructor to create the matrix operation object.
        ///
        /// \param mat_ An **Eigen** sparse matrix object, whose type is
        /// `Eigen::SparseMatrix<Scalar, ...>`.
        ///
        SparseRegularInversePardiso(SparseMatrix& mat_) :
        m_n(mat_.rows()), m_mat(mat_)
        {
            if(mat_.rows() != mat_.cols())
                throw std::invalid_argument("SparseRegularInverse: matrix must be square");
            
            m_pardiso.symbolicFactorization(m_mat);
            m_pardiso.numericalFactorization();
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
        /// Perform the solving operation \f$y=B^{-1}x\f$.
        ///
        /// \param x_in  Pointer to the \f$x\f$ vector.
        /// \param y_out Pointer to the \f$y\f$ vector.
        ///
        // y_out = inv(B) * x_in
        void solve(Scalar* x_in, Scalar* y_out) const
        {
//            MapConstVec x(x_in,  m_n);
//            MapVec      y(y_out, m_n);
//            y.noalias() = m_cg.solve(x);
            
            MapVec x(x_in,  m_n);
            MapVec y(y_out, m_n);
            
            m_pardiso.solve(x);
            y.noalias() = m_pardiso.getX();
        }
        
        ///
        /// Perform the matrix-vector multiplication operation \f$y=Bx\f$.
        ///
        /// \param x_in  Pointer to the \f$x\f$ vector.
        /// \param y_out Pointer to the \f$y\f$ vector.
        ///
        // y_out = B * x_in
        void mat_prod(const Scalar* x_in, Scalar* y_out) const
        {
            MapConstVec x(x_in,  m_n);
            MapVec      y(y_out, m_n);
            y.noalias() = m_mat.template selfadjointView<Uplo>() * x;
        }
    };
    
    
} // namespace Spectra


namespace Gauss {
    
    //Modal Analysis using Spectra
    
    //Temp test Spectre
#include <GenEigsComplexShiftSolver.h>
#include <SymGEigsSolver.h>
#include <GenEigsRealShiftSolver.h>
#include <MatOp/SparseGenMatProd.h>
#include <MatOp/SparseCholesky.h>
#include <MatOp/SparseRegularInverse.h>
#include <MatOp/SparseSymShiftSolve.h>
#include <SymGEigsSolver.h>
#include <stdexcept>
    
    
    
    template<typename DataType, int Flags, typename Indices>
    auto generalizedEigenvalueProblemPardiso(const Eigen::SparseMatrix<DataType, Flags, Indices> &A,
                                                   Eigen::SparseMatrix<DataType, Flags,Indices> &B,
                                                   unsigned int numVecs) {
        
        //Spectra seems to freak out if you use row storage, this copy just ensures everything is setup the way the solver likes
        Eigen::SparseMatrix<DataType> K = -A;
        Eigen::SparseMatrix<DataType, Flags,Indices> M = B;
        
        
        Spectra::SparseSymMatProd<DataType> Aop(K);
        ::Spectra::SparseRegularInversePardiso<DataType>   Bop(M);
        
        //Spectra::SparseSymShiftSolve<DataType> Aop(K);
        
        //Aop.set_shift(1e-3);
        
        // Construct eigen solver object, requesting the smallest three eigenvalues
        Spectra::SymGEigsSolver<DataType, Spectra::SMALLEST_MAGN, Spectra::SparseSymMatProd<DataType>, ::Spectra::SparseRegularInversePardiso<DataType>, Spectra::GEIGS_REGULAR_INVERSE > eigs(&Aop, &Bop, numVecs, 5*numVecs);
        
        
        // Initialize and compute
        eigs.init();
        eigs.compute();
        
        Eigen::VectorXx<DataType> eigsCorrected;
        Eigen::MatrixXx<DataType> evsCorrected; //magnitude of eigenvectors can be wrong in this formulation
        eigsCorrected.resize(eigs.eigenvalues().rows());
        evsCorrected.resize(eigs.eigenvectors().rows(), eigs.eigenvectors().cols());
        
        //int nconv = eigs.compute();
        
        // Retrieve results
        if(eigs.info() == Spectra::SUCCESSFUL) {
            //correct eigenvalues
            for(unsigned int ii=0; ii<eigs.eigenvalues().rows(); ++ii) {
                eigsCorrected[ii] = -(static_cast<DataType>(1)/(eigs.eigenvalues()[ii]));
                evsCorrected.col(ii)  = eigs.eigenvectors().col(ii)/sqrt(eigs.eigenvectors().col(ii).transpose()*M*eigs.eigenvectors().col(ii));
            }
            
            return std::make_pair(eigs.eigenvectors(), eigs.eigenvalues());
        } else {
            std::cout<<"Failure: "<<eigs.info()<<"\n";
            exit(1);
            return std::make_pair(eigs.eigenvectors(), eigs.eigenvalues());
        }
        
    }
    
    template<typename DataType, int Flags, typename Indices>
    auto generalizedEigenvalueProblemPardiso(const Eigen::SparseMatrix<DataType, Flags, Indices> &A,
                                             Eigen::SparseMatrix<DataType, Flags,Indices> &B,
                                             unsigned int numVecs, DataType shift) {
        
        //Spectra seems to freak out if you use row storage, this copy just ensures everything is setup the way the solver likes
        Eigen::SparseMatrix<DataType> K = -A + shift*B;
        Eigen::SparseMatrix<DataType, Flags,Indices> M = B;
        
        
        Spectra::SparseSymMatProd<DataType> Aop(K);
        ::Spectra::SparseRegularInversePardiso<DataType>   Bop(M);
        
        //Spectra::SparseSymShiftSolve<DataType> Aop(K);
        
        //Aop.set_shift(1e-3);
        
        // Construct eigen solver object, requesting the smallest three eigenvalues
        Spectra::SymGEigsSolver<DataType, Spectra::SMALLEST_MAGN, Spectra::SparseSymMatProd<DataType>, ::Spectra::SparseRegularInversePardiso<DataType>, Spectra::GEIGS_REGULAR_INVERSE > eigs(&Aop, &Bop, numVecs, 5*numVecs);
        
        
        // Initialize and compute
        eigs.init();
        eigs.compute();
        
        Eigen::VectorXx<DataType> eigsCorrected;
        Eigen::MatrixXx<DataType> evsCorrected; //magnitude of eigenvectors can be wrong in this formulation
        eigsCorrected.resize(eigs.eigenvalues().rows());
        evsCorrected.resize(eigs.eigenvectors().rows(), eigs.eigenvectors().cols());
        
        //int nconv = eigs.compute();
        
        // Retrieve results
        if(eigs.info() == Spectra::SUCCESSFUL) {
            //correct eigenvalues
            for(unsigned int ii=0; ii<eigs.eigenvalues().rows(); ++ii) {
                eigsCorrected[ii] = -(static_cast<DataType>(1)/(eigs.eigenvalues()[ii]) + shift);
                evsCorrected.col(ii)  = eigs.eigenvectors().col(ii)/sqrt(eigs.eigenvectors().col(ii).transpose()*M*eigs.eigenvectors().col(ii));
            }
            
            return std::make_pair(eigs.eigenvectors(), eigs.eigenvalues());
        } else {
            std::cout<<"Failure: "<<eigs.info()<<"\n";
            exit(1);
            return std::make_pair(eigs.eigenvectors(), eigs.eigenvalues());
        }
        
    }
    
    
}

#endif // SPARSE_REGULAR_INVERSE_PARDISO_H
