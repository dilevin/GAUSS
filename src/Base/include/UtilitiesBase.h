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

template<typename World>
double getEnergy(World &world) {
    
    double energy = 0.0;
    forEach(world.getSystemList(), [&energy, &world](auto a) {
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

//add in the constraints
template<typename Matrix, typename World>
void getConstraintMatrix(Matrix &constraintMatrix, World &world) {
    ASSEMBLEMATINIT(constraintMatrix, world.getNumConstraints(), world.getNumQDotDOFs());
    ASSEMBLELIST(constraintMatrix, world.getConstraintList(), getGradient);
    ASSEMBLEEND(constraintMatrix);
}

//given a a function templated on system type, run it on a system given a system index
struct SystemIndex {
    
    inline SystemIndex() { m_systemType = 0; m_systemIndex = 0; }
    
    inline unsigned int type() const { return m_systemType; }
    inline unsigned int index() const { return m_systemIndex; }
    
    inline unsigned int & index() { return m_systemIndex; }
    inline unsigned int & type() { return m_systemType; }
    
    unsigned int m_systemType;
    unsigned int m_systemIndex;
};



class PassSystem {
    
public: 
    template<typename Func, typename TupleParam, typename ...Params>
    inline auto operator()(TupleParam &tuple, Func &func, SystemIndex &index, Params ...params) {
        return func(tuple[index.index()], params...);
    }
    
};
template<typename SystemList, typename Func, typename ...Params>
inline auto apply(SystemList &list, SystemIndex index, Func &func, Params ...params) {
    PassSystem A;
    apply(list.getStorage(), index.type(), A, func, index, params...);
}


#endif /* UtilitiesBase_h */
