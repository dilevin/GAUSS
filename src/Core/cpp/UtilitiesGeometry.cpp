//
//  UtilitiesGeometry.c
//  Gauss
//
//  Created by David Levin on 5/19/17.
//
//

#include <stdio.h>
#include <UtilitiesEigen.h>
#include <UtilitiesGeometry.h>

//Build Hexahedral element list from array containing regular grid
void Gauss::elementsFromGrid(Eigen::RowVector3i res, Eigen::MatrixXi &F) {
    
    int nX, nY, nZ;

    nX = res(0) - 1;
    nY = res(1) - 1;
    nZ = res(2) - 1;
    
    F.resize(nX*nY*nZ, 8); //resize to build hexahedral elements
    
    auto boxIndex = [&res](unsigned int x, unsigned int y, unsigned int z) { return x + y*res(0)+z*res(0)*res(1); };
    
    //Coords from igl move in x, then y then z
    for(unsigned int ielZ=0; ielZ < nZ; ++ielZ) {
        for(unsigned int ielY=0; ielY < nY; ++ielY) {
            for(unsigned int ielX=0; ielX < nX; ++ielX) {
            
                //element indices are the -x, -y, -z corner of the voxel
                //fill in faces using standard FEM Hex Element ordering
                F.row(ielZ*(nX*nY)+ielY*nX + ielX) <<   boxIndex(ielX, ielY, ielZ),
                                                        boxIndex(ielX+1, ielY, ielZ),
                                                        boxIndex(ielX+1, ielY, ielZ+1),
                                                        boxIndex(ielX, ielY, ielZ+1),
                                                        boxIndex(ielX, ielY+1, ielZ),
                                                        boxIndex(ielX+1, ielY+1, ielZ),
                                                        boxIndex(ielX+1, ielY+1, ielZ+1),
                                                        boxIndex(ielX, ielY+1, ielZ+1);
                
            }
        }
    }
    
}
