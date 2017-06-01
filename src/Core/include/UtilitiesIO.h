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
#include <UtilitiesGeometry.h>

namespace Gauss {

    int openIfstream(std::ifstream &in, std::string filename);
    int openOfstream(std::ofstream &out, std::string filename);
    
    //read in a tetgen file
    void readTetgen(Eigen::MatrixXd &V, Eigen::MatrixXi &F, const std::string nodeFile, const std::string eleFile);
    int loadTet(Eigen::MatrixXd &V, Eigen::MatrixXi &F, std::istream & nodeIn, std::istream & eleIn);
    
    //write data to tetgen
    inline void writeTetgenFiles(const std::string file, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F) { }
    
    //read simulation file (at some point use templates to automatically select the right method, for now these are FEM specific)
    //these are linear hex specific for now (I know, very lazy)
    template<typename Simulation>
    Simulation * readSimFromFile(const std::string filename) {
        
        std::ifstream iFile;
        
        if(openIfstream(iFile, filename) < 0) {
            std::cout<<"Failed to open "<<filename<<"\n";
            return nullptr;
        }
        
        Eigen::MatrixXd V;
        Eigen::MatrixXi F;
        
        
        //read in finite element meshes
        /*From when I was using tetrahedral finite elements
         std::string nodeName, eleName;
        
        iFile >> nodeName;
        iFile >> eleName;
        
        std::cout<<"Reading geometrry from "<<nodeName<<" and "<<eleName<<"\n";
        
        Eigen::MatrixXd V;
        Eigen::MatrixXi F;
        
        readTetgen(V,F, nodeName, eleName);*/
        
        //Hex mesh grids only for now
        //Voxel grid from libigl
        igl::grid(Eigen::RowVector3i(15, 5,5),  V);
        elementsFromGrid(Eigen::RowVector3i(15, 5,5), V, F);

        unsigned int numNodes = 0;
        iFile >> numNodes;
        iFile >> numNodes;
        
        auto *sim = new Simulation(V,F);
        
        //read in and set Young's modulus and poisson's ratio
        double YM, mu;
        YM = 0.0;
        mu = 0.0;
        
        for(unsigned int iel=0; iel<sim->getImpl().getF().rows(); ++iel) {
            iFile>> YM;
            iFile>> mu;
            
            std::cout<<"YM: "<<YM<<" mu: "<<mu<<"\n";
            
            sim->getImpl().getElement(iel)->setParameters(YM, mu);
            sim->getImpl().getElement(iel)->setParameters(YM, mu);
            
            //oFile << sim->getImpl().getElement(iel)->getE() <<" "<<sim->getImpl().getElement(iel)->getMu() <<"\n";
        }

        
        return sim;
        
    }
    
    //write simulation file (see above about specializations)
    template<typename Simulation>
    void writeSimToFile(const std::string simName, const std::string nodeFile, const std::string eleFile, const Simulation *sim) {
        
        //my simple sim setup file format (which will no doubt get more complex as it supports more stuff beyond linear tets
        std::ofstream oFile;
        
        if(openOfstream(oFile, simName+".sim") < 0) {
            return;
        }
        
        //write out geometry files in lines 1 and 2
        //oFile << nodeFile << "\n";
        //oFile << eleFile << "\n";
        
        //total number of elements
        oFile << sim->getImpl().getF().rows()<<"\n";
        oFile << sim->getImpl().getF().cols()<<"\n";
        
        //write out list of element YM and MU
        for(unsigned int iel=0; iel<sim->getImpl().getF().rows(); ++iel) {
            oFile << sim->getImpl().getElement(iel)->getE() <<" "<<sim->getImpl().getElement(iel)->getMu() <<"\n";
        }
        
    }
    
    //read in simulation state
    template<typename DataType>
    void readStateFile(State<DataType> &state, const std::string filename) {
    
        std::ifstream iFile;
        
        if(openIfstream(iFile, filename) < 0) {
            return;
        }
        
        unsigned int stateSize ;
        
        iFile >> stateSize;
        
        std::cout<<"Loading state of size "<<stateSize<<"\n";
    
        if(state.getNumScalarDOF() != stateSize) {
            std::cout<<"Error reading state from disk, state size mismatch \n";
            return;
        }
        
        for(unsigned int ii=0; ii<stateSize; ++ii) {
            iFile>>state[ii];
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
