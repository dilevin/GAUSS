//
//  UtilitiesIO.h
//  Gauss
//
//  Created by David Levin on 4/23/17.
//
//

#ifndef UtilitiesIO_h
#define UtilitiesIO_h

#include <iostream>
#include <fstream>
#include <UtilitiesEigen.h>

namespace Gauss {

    int openIfstream(std::ifstream &in, std::string filename);
    
    //read in a tetgen file
    void readTetgen(Eigen::MatrixXd &V, Eigen::MatrixXi &F, const std::string nodeFile, const std::string eleFile);
    int loadTet(Eigen::MatrixXd &V, Eigen::MatrixXi &F, std::istream & nodeIn, std::istream & eleIn);
}

#endif /* UtilitiesIO_h */
