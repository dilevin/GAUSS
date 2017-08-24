//
//  UtilitiesFEM.h
//  Gauss
//
//  Created by David Levin on 8/23/17.
//
//

#ifndef UtilitiesFEM_h
#define UtilitiesFEM_h

#include <GaussIncludes.h>
#include <FEMIncludes.h>
#include <UtilitiesEigen.h>
#include <complex>

//A collection of useful FEM utilities

namespace Gauss {
    namespace FEM {
        //linear modal analysis
        //returns a pair wherein the first value is a matrix of eigenvectors and the second is vector of corresponding vibrational frequencies
        template<typename World>
        auto linearModalAnalysis(World &world, unsigned int numModes) {
            
            //build mass and stiffness matrices
            AssemblerParallel<double, AssemblerEigenSparseMatrix<double> > mass;
            AssemblerParallel<double, AssemblerEigenSparseMatrix<double> > stiffness;
            
            getMassMatrix(mass, world);
            getStiffnessMatrix(stiffness, world);
            
            auto eigs = generalizedEigenvalueProblem((*stiffness), (*mass), numModes);
            
            //convert to vibrational frequencies
            //hack because generalised eignvalue solver only returns real values
            for(unsigned int ii=0; ii<eigs.second.rows(); ++ii) {
                eigs.second[ii] = std::sqrt(std::fabs(eigs.second[ii]));
            }
            
            return eigs;
        }
    }
}

#endif /* UtilitiesFEM_h */
