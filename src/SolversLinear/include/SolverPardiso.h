//
//  SolverPardiso.h
//  Gauss
//
//  Created by David Levin on 6/5/17.
//
//

#ifndef SolverPardiso_h
#define SolverPardiso_h

#ifdef GAUSS_PARDISO

extern "C" void pardisoinit (void   *, int    *,   int *, int *, double *, int *);
extern "C" void pardiso     (void   *, int    *,   int *, int *,    int *, int *,
 double *, int    *,    int *, int *,   int *, int *,
 int *, double *, double *, int *, double *);
extern "C" void pardiso_chkmatrix  (int *, int *, double *, int *, int *, int *);
extern "C" void pardiso_chkvec     (int *, int *, double *, int *);
extern "C" void pardiso_printstats (int *, int *, double *, int *, int *, int *,
 double *, int *);


#include <cstdio>
#include <cstdlib>
#include <cmath>

template<typename MatrixType>
class SolverPardiso
{
public:
protected:
private:
};

template <>
class SolverPardiso<Eigen::SparseMatrix<double, Eigen::RowMajor> >
{
public:
    
    using SparseMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor>;
    
    //Pardiso stuff
    SolverPardiso(unsigned int numProcessors = 8) {
        
        m_matrixType = 11; //Real symmetric matrix
        
        //Annoying license file stuff
        int error = 0;
        int solver = 0; /* use sparse direct solver */
        pardisoinit (m_ptr,  &m_matrixType, &solver, m_integerParams, m_doubleParams, &error);
        
        if (error != 0)
        {
            if (error == -10 )
                printf("No license file found \n");
            if (error == -11 )
                printf("License is expired \n");
            if (error == -12 )
                printf("Wrong username or hostname \n");
            exit(0);
        }
        else
            printf("[PARDISO]: License check was successful ... \n");
        
        //Number of Processors
        m_integerParams[2] = numProcessors;
        m_maxfct = 1; //perform at most one factorization
        m_mnum = 1; // use factorization 1
        m_nrhs = 1; //1 right hand side for now
        m_displayStats = 0;
        
    }
    
    int symbolicFactorization(Eigen::SparseMatrix<double, Eigen::RowMajor> &A) {
        
        unsigned int nnz = A.nonZeros();
        unsigned int nno = A.outerSize();
        
        m_innerArray.clear();
        m_a.clear();
        for(unsigned int jj=0; jj<nnz; ++jj) {
            m_innerArray.push_back(A.innerIndexPtr()[jj]+1);
            m_a.push_back(A.valuePtr()[jj]);
        }
        
        
        m_outerArray.clear();
        //need to rebuild outer array because Eigen is missing the last index
        for(unsigned int ii=0; ii<nno; ++ii) {
            m_outerArray.push_back(A.outerIndexPtr()[ii] + 1);
        }
        
        m_outerArray.push_back(nnz+1);
        
        //update matrix indices so that they're 1 indexed for the fortran code in pardiso
        int phase = 11;
        double   ddum = 0.;              /* Double dummy */
        int      idum = 0;              /* Integer dummy. */
        int error = 0;
        m_n = A.rows();
        
        /*pardiso_chkmatrix(&m_matrixType, &n, m_a.data(), m_outerArray.data(), m_innerArray.data(), &error);
        if (error != 0) {
            printf("\nERROR in consistency of matrix: %d \n", error);
            exit(1);
        }*/

        pardiso(m_ptr, &m_maxfct, &m_mnum, &m_matrixType, &phase, &m_n, m_a.data(), m_outerArray.data(),
                m_innerArray.data(), &idum, &m_nrhs, m_integerParams, &m_displayStats, &ddum, &ddum, &error, m_doubleParams);
        
        if(error != 0) {
            std::cout<<"Symbolic Factorization Failed: "<<error<<"\n";
            return 0;
        }
        
        return 1;
    }
    
    int numericalFactorization() {
    
        int phase = 22;
        double   ddum = 0.;              /* Double dummy */
        int      idum = 0;              /* Integer dummy. */
        int error = 0;
        
        pardiso(m_ptr, &m_maxfct, &m_mnum, &m_matrixType, &phase, &m_n, m_a.data(), m_outerArray.data(),
                m_innerArray.data(), &idum, &m_nrhs, m_integerParams, &m_displayStats, &ddum, &ddum, &error, m_doubleParams);
        
        if (error != 0) {
            printf("\nERROR during numerical factorization: %d", error);
            return 0;
        }

        return 1;
    }

    template<typename Vector>
    int solve(Vector &rhs) {
        int phase = 33;
        int      idum = 0;              /* Integer dummy. */
        int error = 0;
        
        m_x.resize(rhs.rows(), 1);
        m_integerParams[7] = 1;       /* Max numbers of iterative refinement steps. */
        
        pardiso(m_ptr, &m_maxfct, &m_mnum, &m_matrixType, &phase, &m_n, m_a.data(), m_outerArray.data(),
                m_innerArray.data(), &idum, &m_nrhs, m_integerParams, &m_displayStats, rhs.data(), m_x.data(), &error, m_doubleParams);
        
        /*pardiso (pt, &maxfct, &mnum, &mtype, &phase,
                 &n, a, ia, ja, &idum, &nrhs,
                 iparm, &msglvl, b, x, &error,  dparm); <--- replace dummy variables*/
        
        if (error != 0) {
            printf("\nERROR during solution: %d", error);
            return 0;
        }

        return 1;
    }
    
    template<typename Vector>
    Vector & solve(SparseMatrix &A, Vector &b) {
        symbolicFactorization(A);
        numericalFactorization();
        solve(b);
        cleanup();
        return m_x;
    }
    
    int cleanup() {
        int phase = -1;
        double   ddum = 0.;              /* Double dummy */
        int      idum = 0;              /* Integer dummy. */
        int error = 0;
        
        pardiso(m_ptr, &m_maxfct, &m_mnum, &m_matrixType, &phase, &m_n, m_a.data(), m_outerArray.data(),
                m_innerArray.data(), &idum, &m_nrhs, m_integerParams, &m_displayStats, &ddum, &ddum, &error, m_doubleParams);
        
        assert(error ==0);
        
        return 1;
    }
    
    
    inline auto getX() { return m_x; }
    
protected:
    
    int m_matrixType;
    
    /* Internal solver memory pointer pt,                  */
    /* 32-bit: int pt[64]; 64-bit: long int pt[64]         */
    /* or void *pt[64] should be OK on both architectures  */
    void    *m_ptr[64];
    
    /* Pardiso control parameters. */
    int m_n;
    int      m_integerParams[64];
    double   m_doubleParams[64];
    int      m_maxfct;
    int      m_mnum;
    int      m_nrhs;
    int      m_displayStats;
    std::vector<int> m_outerArray, m_innerArray;
    std::vector<double> m_a;
    Eigen::VectorXd m_x;
    
private:
    
};

#endif

#endif /* SolverPardiso_h */
