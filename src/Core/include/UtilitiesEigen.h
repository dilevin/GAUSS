//
//  UtilitiesEigen.h
//  Gauss
//
//  Created by David Levin on 2/11/17.
//
//

#ifndef UtilitiesEigen_h
#define UtilitiesEigen_h

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Utilities.h>
#include <World.h>


//some useful types
namespace Eigen {
    template<typename DataType>
    using Vector3x = Eigen::Matrix<DataType, 3,1>;
    
    template<typename DataType>
    using VectorXx = Eigen::Matrix<DataType, Eigen::Dynamic, 1>;
    
    template<typename DataType>
    using MatrixXx = Eigen::Matrix<DataType, Eigen::Dynamic, 1>;
}

namespace Gauss {
    //state ptr direct to eigen map for a single property (position or velocity)
    template<unsigned int Property, typename World>
    Eigen::Map<Eigen::VectorXd> mapStateEigen(World &world) {
        std::pair<double *, unsigned int> ptr = world.getState().template getStatePtr<Property>();
        return Eigen::Map<Eigen::VectorXd>(ptr.first, ptr.second);
    }

    //state ptr for the whole thing
    template<typename World>
    Eigen::Map<Eigen::VectorXd> mapStateEigen(World &world) {
        std::pair<double *, unsigned int> ptr = world.getState().getStatePtr();
        return Eigen::Map<Eigen::VectorXd>(ptr.first, ptr.second);
    }
    
    template<typename DOF, typename DataType, typename ...SystemTypes, typename ...ForceTypes, typename ...ConstraintTypes>
    inline Eigen::Map<Eigen::VectorXd> mapDOFEigen(DOF &dof, World<DataType,
                                                   std::tuple<SystemTypes...>,
                                                   std::tuple<ForceTypes...>,
                                                   std::tuple<ConstraintTypes...> > &world) {
        std::pair<double *, unsigned int> qPtr = dof.getPtr(world.getState());
        //set position DOF and check
        return Eigen::Map<Eigen::VectorXd>(qPtr.first, dof.getNumScalarDOF());

    }
    
    template<typename DOF, typename DataType>
    inline Eigen::Map<Eigen::VectorXd> mapDOFEigen(DOF &dof, const State<DataType> &state) {
        std::pair<double *, unsigned int> qPtr = dof.getPtr(state);
        //set position DOF and check
        return Eigen::Map<Eigen::VectorXd>(qPtr.first, dof.getNumScalarDOF());
        
    }

}
#endif /* UtilitiesEigen_h */
