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
#include <DOFParticle.h>
#include <DOFRotation.h>
#include <DOFPair.h>
#include <DOFList.h>
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

template<typename Matrix, typename World>
void getInternalForceVector(Matrix &forceVector, World &world) {
    ASSEMBLEVECINIT(forceVector, world.getNumQDotDOFs());
    ASSEMBLELIST(forceVector, world.getSystemList(), getInternalForce);
    ASSEMBLEEND(forceVector);
}


//get strain energy
template<typename World>
 double getStrainEnergy(World &world) {
    
     double energy = 0.0;
     
     forEach(world.getSystemList(), [&energy, &world](auto a) {
         energy += a->getStrainEnergy(world.getState());
     });
     
     return energy;
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

//Specific map for rotations
template<typename DataType, unsigned int Property>
inline Eigen::Map<Eigen::Quaternion<DataType> > mapDOFEigenQuat(const DOFRotation<DataType, Property> &dof, const State<DataType> &state) {
    std::tuple<double *, unsigned int> qPtr = dof.getPtr(state);
    return Eigen::Map<Eigen::Quaternion<DataType> >(std::get<0>(qPtr));
}

//Initializers for DOFS
//Default Initializers just zeros things out
template<typename DOFType>
class InitializeDOFClass
{
    
public:
    template<typename State>
    explicit inline InitializeDOFClass(DOFType &dof, State &state) {
        std::cout<<"Should not be here \n";
        exit(0);
    }
};

template<typename DataType, unsigned int PropertyIndex>
class InitializeDOFClass<DOFRotation<DataType,PropertyIndex> >
{
    
public:
    template<typename State>
    explicit inline InitializeDOFClass(DOFRotation<DataType,PropertyIndex> &dof, State &state) {
        std::cout<<"Rotation \n";
        
        auto statePtr = dof.getPtr(state);
        std::memset(std::get<0>(statePtr), 0, sizeof(DataType)*std::get<1>(statePtr));
        std::get<0>(statePtr)[3] = 1.0;
    }
};

template<typename DataType, unsigned int PropertyIndex>
class InitializeDOFClass<DOFParticle<DataType,PropertyIndex> >
{
    
public:
    template<typename State>
    explicit inline InitializeDOFClass(DOFParticle<DataType,PropertyIndex> &dof, State &state) {
        std::cout<<"Position \n";
        //standard initializer sets everything to zero
        auto statePtr = dof.getPtr(state);
        std::memset(std::get<0>(statePtr), 0, sizeof(DataType)*std::get<1>(statePtr));
    }
};

template<typename DataType, unsigned int PropertyIndex, template<typename A, unsigned int B> class DOF1, template<typename A, unsigned int B> class DOF2>
class InitializeDOFClass< DOFPair<DataType, DOF1, DOF2, PropertyIndex> >
{
    
public:
    template<typename State>
    explicit inline InitializeDOFClass(DOFPair<DataType,DOF1, DOF2, PropertyIndex> &dof, State &state) {
        std::cout<<"DOFPair \n";
        InitializeDOFClass<DOF1<DataType, PropertyIndex> >(dof.first(), state);
        InitializeDOFClass<DOF2<DataType, PropertyIndex>>(dof.second(), state);
    }
};

template<typename DOF, typename DataType>
inline void InitializeDOF(DOF &dof, State<DataType> &state) {
    InitializeDOFClass<DOF>(dof, state);
}

//Initialize everything
template<typename World>
void initializeDOFs(World &world) {
    forEach(world.getSystemList(), [&world](auto a){
        InitializeDOF(a->getQ(), world.getState());
        InitializeDOF(a->getQDot(), world.getState());
    });
}




//incrementing DOFs
template<typename QDOF, typename QDOTDOF, typename DataType>
class IncrementDOFClass
{
    
public:
    template<typename State>
    explicit inline IncrementDOFClass(QDOF &q, QDOTDOF &qDot, DataType a, State &state) {
        
        //normal addition
        #pragma omp parallel for
        for(unsigned int ii=0; ii<q.getNumScalarDOF(); ++ii) {
            std::get<0>(q.getPtr(state))[ii] += a*std::get<0>(qDot.getPtr(state))[ii];
        }
    }
};

//deal with rotations (standard addition doesn't work)
template<unsigned int PropertyIndex, typename DataType>
class IncrementDOFClass<DOFRotation<DataType,PropertyIndex>, DOFParticle<DataType,PropertyIndex>, DataType >
{
    
public:
    template<typename State>
    explicit inline IncrementDOFClass(DOFRotation<DataType,PropertyIndex> &q, DOFParticle<DataType,PropertyIndex> &qDot, DataType a, State &state) {
        
        std::cout<<"Update Rotation \n";
        
        //convert angular velocity to quaternion and post multiply to update current rotation
        mapEigenDOFEigenQuat(q, state) = Eigen::Quaternion<DataType>(Eigen::AngleAxis<DataType>(a*mapEigenDOF(qDot, state).norm(), mapEigenDOF(qDot, state).normalized()))*mapEigenDOFEigenQuat(q, state);
    }
};

template<typename QDOF, typename QDOTDOF, typename DataType, typename State>
inline void incrementDOF(QDOF &q, QDOTDOF &qDot, DataType a, State &state) {
    IncrementDOFClass<QDOF, QDOTDOF, DataType>(q, qDot, a, state);
}

//Update is a covenience method to do the following operation that occurs all the time
// q = q + dt*qDot where + is approriate to the particular DOF
template<typename World, typename State, typename DataType>
inline void updateState(World &world, State & state, DataType dt) {
    
    std::cout<<"UPDATE STATE\n";
    //update position in state
    forEach(world.getSystemList(), [&world, &state, &dt](auto a) {
       
        #pragma omp task shared(world, state)
        {
            //iterate through dofs in this system and do the update
            incrementDOF(a->getQ(), a->getQDot(), dt, state);
        }
        
    });
    
    std::cout<<"DONE UPDATE STATE\n";
    
}
#endif /* UtilitiesBase_h */
