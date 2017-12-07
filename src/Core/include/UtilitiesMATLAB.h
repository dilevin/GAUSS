//
//  UtilitiesMATLAB.h
//  Gauss
//
//  Created by David Levin on 4/17/17.
//
//

#include <iostream>
#include <string>
#include <UtilitiesEigen.h>
#include <unsupported/Eigen/SparseExtra>

//Code that is helpful for interfacing with MATLAB
#ifndef UtilitiesMATLAB_h
#define UtilitiesMATLAB_h

namespace Gauss {
    
    //Write matrix types for matlab
    template<typename MatrixType>
    inline void toMatlab(MatrixType &toWrite, std::string matrixName) {
        // Init
        std::stringstream stream;
        std::cout<<"Write to matlab, generic method == something has gone horribly wrong :( \n";
        /*
         list<pair<size_t, std::complex<double> > >::iterator it;
         list<pair<size_t, std::complex<double> > >::iterator end;
         
         // Common part
         stream << matlabCommon(matrixName);
         
         // Values
         stream << "[";
         for(size_t i = 0; i < nRow; i++){
         it  = data[i].begin();
         end = data[i].end();
         
         for(; it != end; it++)
         stream << std::scientific << std::setprecision(16)
         << "(" << it->second.real() << " + i * "<< it->second.imag() << ")"
         << ", ";
         }
         stream << "], ";
         
         // Number of rows and columns
         stream << nRow << ", " << nCol << ");";
         */
        // Return
        
    }

    template<>
    inline void toMatlab<Eigen::SparseMatrix<double> >(Eigen::SparseMatrix<double> &toWrite, std::string filename) {
        std::cout<<"Writing Eigen Sparse Matrix to MATLAB \n";
        
        Eigen::saveMarket(toWrite, filename);

        
        
    }

    template<>
    inline void toMatlab<Eigen::SparseMatrix<double, Eigen::RowMajor> >(Eigen::SparseMatrix<double,Eigen::RowMajor> &toWrite, std::string filename) {
        std::cout<<"Writing Eigen Sparse Matrix to MATLAB \n";
        
        Eigen::saveMarket(toWrite, filename);
        
        
        
    }

    
    template<>
    inline void toMatlab<Eigen::Matrix<double, 12, 12> >(Eigen::Matrix<double, 12,12> &toWrite, std::string filename) {
        std::ofstream file(filename);
        if (file.is_open())
        {
            file << toWrite;
        }
        file.close();
    }
    
    template<>
    inline void toMatlab<Eigen::MatrixXd>(Eigen::MatrixXd &toWrite, std::string filename) {
        std::ofstream file(filename);
        if (file.is_open())
        {
            file << toWrite;
        }
        file.close();
    }

    
    template<>
    inline void toMatlab<Eigen::VectorXd>(Eigen::VectorXd &toWrite, std::string filename) {
        std::ofstream file(filename);
        if (file.is_open())
        {
            file << toWrite;
        }
        file.close();
    }
    
    template<>
    inline void toMatlab<Eigen::Map<Eigen::VectorXd> >(Eigen::Map<Eigen::VectorXd> &toWrite, std::string filename) {
        std::ofstream file(filename);
        if (file.is_open())
        {
            file << toWrite;
        }
        file.close();
    }
}

#endif /* UtilitiesMATLAB_h */
