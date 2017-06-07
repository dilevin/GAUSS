
#ifndef _WORLD_H
#define _WORLD_H

#include "CoreDefines.h"
#include "CoreOSSpecific.h"

#include "State.h"
#include "MultiVector.h"
#include <tuple>

/*
 * The world holds the lists of all systems as well as the list of all forces and constraints that act on those systems
 * The world functions as an iterator for these lists
 */
namespace Gauss
{
    template<typename DataType, typename Impl> class PhysicalSystem;
    template<typename Datatype, typename Impl > class Force;
    template<typename Datatype, typename Impl > class Constraint;
    
    template <typename DataType, typename SystemTuple, typename ForceTuple, typename ConstraintTuple> class World;
    
    template<typename DataType, typename ...SystemTypes, typename ...ForceTypes, typename ...ConstraintTypes>
    class World<DataType, std::tuple<SystemTypes...>, std::tuple<ForceTypes...>, std::tuple<ConstraintTypes...> >
    {
        
    public:
        explicit World();
        explicit World(const World &toCopy);
        
        ~World();
        
        //Accessors
        inline double getCurrentTime() const { return m_t; }
        inline void setCurrentTime(const double time) { m_t = time; }
        
        inline State<DataType> & getState() { return m_state; }
        inline const State<DataType> & getState() const { return m_state; }
        
        inline const MultiVector<SystemTypes...> & getSystemList() const { return m_systems; }
        inline MultiVector<SystemTypes...> & getSystemList() { return m_systems; }
        
        inline const MultiVector<ForceTypes...> & getForceList() const { return m_forces; }
        inline MultiVector<ForceTypes...> & getForceList() { return m_forces; }
        
        inline const MultiVector<ConstraintTypes...> & getConstraintList() const { return m_constraints; }
        inline MultiVector<ConstraintTypes...> & getConstraintList() { return m_constraints; }
        
        
        //Add items to the world
        template<typename Impl>
        int addSystem(PhysicalSystem<DataType, Impl> *system);
        
        template<typename Impl>
        int addForce(Force<DataType, Impl> *force);
        
        template<typename Impl>
        int addConstraint(Constraint<DataType, Impl> *constraint);
        
        //Get globaal properties of the worls
        inline unsigned int getNumQDOFs() const { return m_numQDOFs; }
        inline unsigned int getNumQDotDOFs() const { return m_numQDotDOFs; }
        inline unsigned int getTotalDOFs() const { return m_numQDOFs + m_numQDotDOFs; }
        
        //Matrix Storage Stuff
        inline unsigned int getMassNNZ()  {
            unsigned int nnz = 0;
            forEach(m_systems, [&nnz](auto a) {
                nnz += a->getMassNNZ();
            });
            return nnz;
        }
        
        inline unsigned int getStiffnessNNZ()  {
            unsigned int nnz = 0;
            forEach(m_systems, [&nnz](auto a) {
                nnz += a->getStiffnessNNZ();
            });
            
            forEach(m_forces, [&nnz](auto a) {
                nnz += a->getStiffnessNNZ();
            });
            
            return nnz;
        }
        
        inline unsigned int getConstraintNNZ()  {
            unsigned int nnz = 0;
            forEach(m_constraints, [&nnz](auto a) {
                nnz += a->getNNZ();
            });
            return nnz;
        }

        
        
        
        
        
        //constraint stuff
        inline unsigned int getNumConstraints() const { return m_numConstraints; }
        
        //Check memory make sure everything is cool before kicking off calcultions
        int finalize();
        
    protected:
        
        
        double m_t; //world time
        
        //Lists
        MultiVector<SystemTypes...> m_systems;
        MultiVector<ForceTypes...> m_forces;
        MultiVector<ConstraintTypes...> m_constraints;
        
        unsigned int m_numQDOFs;
        unsigned int m_numQDotDOFs;
        unsigned int m_numConstraints;
        
        State<DataType> m_state;
        
    private:
    };
}

#include "World.h"

using namespace Gauss;

template<typename DataType, typename ...SystemTypes, typename ...ForceTypes, typename ...ConstraintTypes>
World<DataType, std::tuple<SystemTypes...>, std::tuple<ForceTypes...>, std::tuple<ConstraintTypes...> >::World() : m_state(0)
{
    m_t = 0.0;
    m_numQDOFs = 0;
    m_numQDotDOFs = 0;
    m_numConstraints = 0;
    m_numConstraints = 0;
}

template<typename DataType, typename ...SystemTypes, typename ...ForceTypes, typename ...ConstraintTypes>
World<DataType,std::tuple<SystemTypes...>, std::tuple<ForceTypes...>, std::tuple<ConstraintTypes...> >::World(const World &/*toCopy*/)
{
    DEBUGPRINT("World: Copy constructor does nothing");
}

template<typename DataType, typename ...SystemTypes, typename ...ForceTypes, typename ...ConstraintTypes>
World<DataType,std::tuple<SystemTypes...>, std::tuple<ForceTypes...>, std::tuple<ConstraintTypes...> >::~World()
{
    DEBUGPRINT("World: Destructor deletes all objects");
    m_systems.clear();
    m_forces.clear();
    //m_constraints.erase(m_constraints.begin(), m_constraints.end());
    
}

template<typename DataType, typename ...SystemTypes, typename ...ForceTypes, typename ...ConstraintTypes>
template<typename Impl>
int World<DataType,std::tuple<SystemTypes...>, std::tuple<ForceTypes...>, std::tuple<ConstraintTypes...> >::addSystem(PhysicalSystem<DataType,Impl> *system)
{
  m_systems.add(system) ;
  return 1;
}

template<typename DataType, typename ...SystemTypes, typename ...ForceTypes, typename ...ConstraintTypes>
template<typename Impl>
int World<DataType,std::tuple<SystemTypes...>, std::tuple<ForceTypes...>, std::tuple<ConstraintTypes...> >::addForce(Force<DataType, Impl> *force)
{
    m_forces.add(force);
    return 1;
}

template<typename DataType, typename ...SystemTypes, typename ...ForceTypes, typename ...ConstraintTypes>
template<typename Impl>
int World<DataType,std::tuple<SystemTypes...>, std::tuple<ForceTypes...>, std::tuple<ConstraintTypes...> >::addConstraint(Constraint<DataType, Impl> *constraint)
{
    m_constraints.add(constraint);
    return 1;
}

template<typename DataType, typename ...SystemTypes, typename ...ForceTypes, typename ...ConstraintTypes>
int World<DataType, std::tuple<SystemTypes...>, std::tuple<ForceTypes...>, std::tuple<ConstraintTypes...> >::finalize()
{
    //start
    unsigned int totalQDotDOFs = 0;
    unsigned int totalQDOFs = 0;
    
    //compute the total number of Dofs in the system
    forEach(m_systems, [&totalQDOFs, &totalQDotDOFs](auto a){totalQDOFs+=a->getQ().getNumScalarDOF();
                                                         totalQDotDOFs+=a->getQDot().getNumScalarDOF();});
    std::cout<<"Number of DOFS: "<<(totalQDOFs+totalQDotDOFs)<<"\n";
    
    //allocate memory in the state
    m_state.resize(totalQDOFs+totalQDotDOFs);
    m_state.setOffset(totalQDOFs);
    m_numQDOFs = totalQDOFs;
    m_numQDotDOFs = totalQDotDOFs;
    //update DOF global indices
    unsigned int currentOffsetQ = 0;
    unsigned int currentOffsetQDot = 0;
    forEach(m_systems, [&currentOffsetQ, &currentOffsetQDot, &totalQDOFs](auto a){a->getQ().offsetGlobalId(currentOffsetQ);
                                                a->getQDot().offsetGlobalId(currentOffsetQDot);
                                                currentOffsetQ += a->getQ().getNumScalarDOF();
                                                currentOffsetQDot += a->getQDot().getNumScalarDOF(); });
    
    //done
    
    
    //update constraint global indices (just count out constraints)
    unsigned int totalConstraints = 0;
    forEach(m_constraints, [&totalConstraints](auto a) {
        a->getIndex().offsetGlobalId(totalConstraints);
        totalConstraints += a->getIndex().getNumScalarDOF();
    });
    
    m_numConstraints = totalConstraints;
    
    return 1;
}

#endif // WORLD_H
