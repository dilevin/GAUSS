//
//  UtilitiesGeometry.h
//  Gauss
//
//  Created by David Levin on 5/19/17.
//
//

#ifndef UtilitiesGeometry_h
#define UtilitiesGeometry_h

#include <UtilitiesEigen.h>
#include <igl/grid.h>

namespace Gauss {

    void elementsFromGrid(Eigen::RowVector3i res, Eigen::MatrixXd &V, Eigen::MatrixXi &F);
    
}

#endif /* UtilitiesGeometry_h */
