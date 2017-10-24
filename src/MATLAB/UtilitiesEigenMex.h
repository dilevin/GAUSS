//
//  UtilitiesEigenMex.h
//  Gauss
//
//  Created by David Levin on 9/11/17.
//
//

#ifndef UtilitiesEigenMex_h
#define UtilitiesEigenMex_h

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "mex.h"
#include "matrix.h"

Eigen::MatrixXd matlabToDouble(const mxArray *matlabMatrix) {
    
    //check array type
    if(!mxIsDouble(matlabMatrix)) {
        mexErrMsgTxt("Matrix is not of type double");
    }
    
    unsigned int m = mxGetM(matlabMatrix);
    unsigned int n = mxGetN(matlabMatrix);
    
    Eigen::MatrixXd data(m,n);
    data.setZero();
    for(unsigned int ii=0; ii<m; ++ii) {
        for(unsigned int jj=0; jj<n;++jj) {
            data(ii,jj) = (mxGetPr(matlabMatrix))[ii+m*jj];
        }
    }
    
    return data;
}

Eigen::MatrixXi matlabToInt32(const mxArray *matlabMatrix) {
    
    //check array type
    if(!mxIsInt32(matlabMatrix)) {
        mexErrMsgTxt("Matrix is not of type Int64");
    }
    
    unsigned int m = mxGetM(matlabMatrix);
    unsigned int n = mxGetN(matlabMatrix);
    
    Eigen::MatrixXi data(m,n);
    data.setZero();
    for(unsigned int ii=0; ii<m; ++ii) {
        for(unsigned int jj=0; jj<n;++jj) {
            data(ii,jj) = ((int *)mxGetData(matlabMatrix))[ii+m*jj]-1;
        }
    }
    
    return data;
}

template<typename DataType, int Options>
mxArray * eigenSparseToMATLAB(const Eigen::SparseMatrix<DataType, Options> &matrix) {
    
    mexPrintf("Base sparse matrix data copy");
    //unsigned int nnz = matrix.nonZeros();
    
}

template<>
mxArray * eigenSparseToMATLAB<double,  Eigen::ColMajor>(const Eigen::SparseMatrix<double, Eigen::ColMajor> &matrix) {
    
    unsigned int nnz = matrix.nonZeros();
    
    mxArray *sparseArray = mxCreateSparse(matrix.rows(), matrix.cols(), nnz, mxREAL);
    
    double *pr = mxGetPr(sparseArray);
    mwIndex *ir=mxGetIr(sparseArray);
    mwIndex *jc = mxGetJc(sparseArray);
    
    //pretty straight conversion, just need to add one to my indices
    for(unsigned int ii=0; ii<nnz; ++ii) {
        //set actual values and row indices
        pr[ii] = matrix.valuePtr()[ii];
        ir[ii] = matrix.innerIndexPtr()[ii];
    }
    
    //set column indices
    for(unsigned int ii=0; ii<matrix.cols()+1;++ii) {
        jc[ii] = matrix.outerIndexPtr()[ii];
    }
        
    return sparseArray;
    
}

//NOTE YOU NEED TO TRANSFER THE MATRIX ON WRITE OUT
template<>
mxArray * eigenSparseToMATLAB<double,  Eigen::RowMajor>(const Eigen::SparseMatrix<double, Eigen::RowMajor> &matrix) {
    
    unsigned int nnz = matrix.nonZeros();
    
    mxArray *sparseArray = mxCreateSparse(matrix.rows(), matrix.cols(), nnz, mxREAL);
    
    double *pr = mxGetPr(sparseArray);
    mwIndex *ir=mxGetIr(sparseArray);
    mwIndex *jc = mxGetJc(sparseArray);
    
    mwIndex *w = new mwIndex[matrix.cols()+1];
    memset(w, 0, sizeof(mwIndex)*(matrix.cols()+1));
    
    for(unsigned int ii=0; ii<nnz; ++ii) {
        w[matrix.innerIndexPtr()[ii]+1]++;
    }
    
    for(unsigned int ii=0; ii<matrix.cols(); ++ii) {
        w[ii+1] = w[ii]+w[ii+1];
        jc[ii+1] = w[ii+1];
        
        //mexPrintf("%zu \n", w[ii+1]);
    }
    
    jc[0] = 0;
    
    mwIndex q = 0;
    
    for(unsigned int jj=0; jj<matrix.rows(); ++jj) {
        
        for(unsigned int pp = matrix.outerIndexPtr()[jj]; pp < matrix.outerIndexPtr()[jj+1]; ++pp) {
            //ir[index] = column where index is the
            q = w[matrix.innerIndexPtr()[pp]]++;
            ir[q] = jj;
            pr[q] = matrix.valuePtr()[pp];
        }
    }
    
    return sparseArray;
    
    
}

template< int Rows,  int Cols>
mxArray * eigenDenseToMATLAB(const Eigen::Matrix<double, Rows,Cols> &matrix) {
    //mexPrintf("copy eigen dense matrix to matlab\n");
    
    mwSize dims[2];
    dims[0] = matrix.rows();
    dims[1] = matrix.cols();
    mxArray *denseArray = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
    memcpy(mxGetPr(denseArray), matrix.data(), sizeof(double)*matrix.rows()*matrix.cols());
    return denseArray;
}


#endif /* UtilitiesEigenMex_h */
