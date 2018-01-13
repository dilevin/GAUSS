//
//  UtilitiesBase.h
//  Gauss
//
//  Created by David Levin on 6/1/17.
//
// Some useful methods for dealing with aggregating data across physical systems and what not
#ifndef UtilitiesBase_h
#define UtilitiesBase_h

#include <Assembler.h>
#include <PhysicalSystem.h>
#include <igl/boundary_facets.h>
#include <igl/writeOBJ.h>

template<typename World>
double getEnergy(World &world) {
    
    double energy = 0.0;
    forEach(world.getSystemList(), [&energy, &world](auto a) {
        energy += a->getEnergy(world.getState());
    });
    
    forEach(world.getForceList(), [&energy, &world](auto a) {
        energy += a->getEnergy(world.getState());
    });
    
    return energy;
}

template<typename Matrix, typename World>
void getMassMatrix(Matrix &massMatrix, World &world) {

    //get mass matrix
    ASSEMBLEMATINIT(massMatrix, world.getNumQDotDOFs(), world.getNumQDotDOFs());
    ASSEMBLELIST(massMatrix, world.getSystemList(), getMassMatrix);
    ASSEMBLEEND(massMatrix);
}

template<typename Matrix, typename World>
void getStiffnessMatrix(Matrix &stiffnessMatrix, World &world) {
 
    //get stiffness matrix
    ASSEMBLEMATINIT(stiffnessMatrix, world.getNumQDotDOFs(), world.getNumQDotDOFs());
    ASSEMBLELIST(stiffnessMatrix, world.getSystemList(), getStiffnessMatrix);
    ASSEMBLELIST(stiffnessMatrix, world.getForceList(), getStiffnessMatrix);
    ASSEMBLEEND(stiffnessMatrix);

    
}

template<typename Matrix, typename World>
void getForceVector(Matrix &forceVector, World &world) {
    ASSEMBLEVECINIT(forceVector, world.getNumQDotDOFs());
    ASSEMBLELIST(forceVector, world.getForceList(), getForce);
    ASSEMBLELIST(forceVector, world.getSystemList(), getForce);
    ASSEMBLEEND(forceVector);
}

template<typename Matrix, typename System, typename World>
void getForceVector(Matrix &forceVector, System &system, World &world) {
    ASSEMBLEVECINIT(forceVector, system.getQ().getNumScalarDOF());
    forceVector.setOffset(-system.getQ().getGlobalId(), 0);
    system.getForce(forceVector, world.getState());
    //ASSEMBLELIST(forceVector, world.getForceList(), getForce);
    //ASSEMBLELIST(forceVector, world.getSystemList(), getForce);
    ASSEMBLEEND(forceVector);
}

template<typename Matrix, typename System, typename World>
void getInternalForceVector(Matrix &forceVector, System &system, World &world) {
    ASSEMBLEVECINIT(forceVector, system.getQ().getNumScalarDOF());
    forceVector.setOffset(-system.getQ().getGlobalId(), 0);
    system.getInternalForce(forceVector, world.getState());
    //ASSEMBLELIST(forceVector, world.getForceList(), getForce);
    //ASSEMBLELIST(forceVector, world.getSystemList(), getForce);
    ASSEMBLEEND(forceVector);
}


//add in the constraints
template<typename Matrix, typename World>
void getConstraintMatrix(Matrix &constraintMatrix, World &world) {
    ASSEMBLEMATINIT(constraintMatrix, world.getNumConstraints(), world.getNumQDotDOFs());
    ASSEMBLELIST(constraintMatrix, world.getConstraintList(), getGradient);
    ASSEMBLEEND(constraintMatrix);
}

//given a a function templated on system type, run it on a system given a system index
struct SystemIndex {
    
    inline SystemIndex() {
        m_type = -1;
        m_index = 0;
    }
    
    inline SystemIndex(unsigned int type, unsigned int index) {
        m_type = type;
        m_index = index;
    }
    
    inline int & type()  { return m_type; }
    inline int & index() { return m_index; }
    
    inline const int & type() const { return m_type; }
    inline const int & index() const { return m_index; }
    
    int m_type; //-1 is fixed object, doesn't need collision response
    int m_index; //index of object in respective systems list
    
    
};


class PassSystem {
    
public:
    template<typename Func, typename TupleParam, typename ...Params>
    inline decltype(auto) operator()(TupleParam &tuple, Func &func, SystemIndex &index, Params ...params) {
        return func(tuple[index.index()], params...);
    }
    
};

template<typename SystemList, typename Func, typename ...Params>
inline decltype(auto) apply(SystemList &list, SystemIndex index, Func &func, Params ...params) {
    PassSystem A;
    apply(list.getStorage(), index.type(), A, func, index, params...);
}

template<typename SystemList, typename Func, typename ...Params>
inline decltype(auto) apply(SystemList &list, SystemIndex index, Func func, Params ...params) {
    PassSystem A;
    apply(list.getStorage(), index.type(), A, func, index, params...);
}

template<typename Geometry>
inline void writeGeoToFile(std::string filename, Geometry &geo, Eigen::VectorXd &u) {
    std::cout<<"This write GEO method does nothing\n";
}

template<>
inline void writeGeoToFile<std::pair<Eigen::MatrixXd &, Eigen::MatrixXi &> >(std::string filename, std::pair<Eigen::MatrixXd &, Eigen::MatrixXi &> &geo, Eigen::VectorXd &u) {
    
    Eigen::MatrixXi B; //boundary facets
    Eigen::MatrixXd  uMat = Eigen::Map<Eigen::MatrixXd>(u.data(), 3, u.rows()/3);
    
    std::cout<<"Writing "<<filename<<"\n";
    
    //get the boundary facets for my data then write everything to disk
    igl::boundary_facets(geo.second, B);
    
    B = B.rowwise().reverse().eval();
    
    igl::writeOBJ(filename, geo.first+uMat.transpose(), B);

}

//write obj file for each object in scene something like 'simname_objindex_frame_index.obj'
template<typename World>
inline void writeWorldToOBJ(std::string folder, std::string simName, World &world, unsigned int frameNumber) {
    
    //iterate through world, get geometry for each system and write to OBJ
    std::cout<<"WARNING Only works for FEM Systems Currently\n";
    
    //build protostring for file names
    std::string firstPart = folder+"/"+simName;
    unsigned int numObjects = world.getNumSystems();
    
    //Loop through every object, check if any points are on the wrong side of the floor, if so
    //record collision
    forEachIndex(world.getSystemList(), [&world, &firstPart, &numObjects, &frameNumber](auto type, auto index, auto &a) {
        
        auto geo = a->getGeometry();
        
        int objID = type*numObjects + index;
        
        std::string padFrameNumber = std::string(10-std::to_string(frameNumber).size(), '0').append(std::to_string(frameNumber));
        
        std::string outputFile = firstPart + "_"+std::to_string(objID)+"_"+padFrameNumber+".obj";
    
        //get object displacuments
        Eigen::VectorXd disp = mapDOFEigen(a->getQ(), world.getState());
        
        writeGeoToFile(outputFile, geo, disp);
    
        
       
    });
    
    
}
#endif /* UtilitiesBase_h */
