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
    int openOfstream(std::ofstream &out, std::string filename);
    
    //read in a tetgen file
    void readTetgen(Eigen::MatrixXd &V, Eigen::MatrixXi &F, const std::string nodeFile, const std::string eleFile);
    int loadTet(Eigen::MatrixXd &V, Eigen::MatrixXi &F, std::istream & nodeIn, std::istream & eleIn);
    
    //write data to tetgen
    inline void writeTetgenFiles(const std::string file, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F) { }
    
    //write simulation file
    template<typename Simulation>
    void writeSimToFile(const std::string simName, const std::string nodeFile, const std::string eleFile, const Simulation *sim) {
        
        //my simple sim setup file format (which will no doubt get more complex as it supports more stuff beyond linear tets
        std::ofstream oFile;
        
        if(openOfstream(oFile, simName+".sim") < 0) {
            return;
        }
        
        //write out geometry files in lines 1 and 2
        oFile << nodeFile << "\n";
        oFile << eleFile << "\n";
        
        //total number of elements
        oFile << sim->getImpl().getF().rows()<<"\n";
        
        //write out list of element YM and MU
        for(unsigned int iel=0; iel<sim->getImpl().getF().rows(); ++iel) {
            oFile << sim->getImpl().getElement(iel)->getE() <<" "<<sim->getImpl().getElement(iel)->getMu() <<"\n";
        }
        
    }
    
    //write out simulation state
    //line 1: numQ, numQDot
    //line 2: through numQ: dofs
    //line numQ+1 .... qdots
    template<typename DataType>
    void writeStateToFile(const std::string simName, State<DataType> &state) {
        
        std::ofstream oFile;
        
        if(openOfstream(oFile, simName+".traj") < 0) {
            return;
        }
        
        oFile<<state.getNumScalarDOF()<<"\n";
        
        for(unsigned int ii=0; ii<state.getNumScalarDOF(); ++ii) {
           oFile<<state[ii]<<"\n";
        }
        
    }
}



#endif /* UtilitiesIO_h */
