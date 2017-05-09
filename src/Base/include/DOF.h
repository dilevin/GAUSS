#ifndef DOF_H
#define DOF_H

#include "State.h"

namespace Gauss {

    /**
     Even baser class
     */
    template<typename DataType, unsigned int Property=0>
    class DOFBase
    {
    public:
        explicit DOFBase(unsigned int localId) {
            m_localId = localId;
            m_globalId = localId;
        }
        
        ~DOFBase() { }
        
        inline void setIds(unsigned int localId, unsigned int globalId) {
            m_localId = localId;
            m_globalId = globalId;
        }
        
        inline unsigned int getNumScalarDOF() {
            return m_numScalarDOF;
        }
        
        inline unsigned int getNumScalarDOF() const {
            return m_numScalarDOF;
        }
        
        inline unsigned int getLocalId() const {
            return m_localId;
        }
        
        inline unsigned int getGlobalId() const {
            return m_globalId;
        }
        
        inline std::tuple<DataType *, unsigned int> getPtr(const State<DataType> &state, unsigned int offset) {
            return state. template getStatePtr<Property>(offset);
        }
        
        inline std::tuple<DataType *, unsigned int> getPtr(const State<DataType> &state) {
            return state. template getStatePtr<Property>(m_globalId);
        }
        
        virtual inline void offsetGlobalId(unsigned int offset) {
            m_globalId += offset;
        }
        
        
    protected:
        
        unsigned int m_localId;
        unsigned int m_globalId;
        unsigned int m_numScalarDOF;
        
    private:
    };
    
    /**
     "Base Class" for Degree-Of-Freedom for a physical system
     DOFs serve as a reference to memory in the state and can contain methods that manipulate that data. They do not store there own information
     */
    template<typename DataType, template<typename A, unsigned int B> class Impl, unsigned int PropertyIndex = 0>
    class DOF : public DOFBase<DataType, PropertyIndex>
    {
    public:
        
        //Do something about this, dofID should be pulled from the state or assigned automatically
        //somehow
        explicit DOF(unsigned int localId = 0) : DOFBase<DataType, PropertyIndex>(localId), m_dofImpl() {
            DOFBase<DataType, PropertyIndex>::m_numScalarDOF = m_dofImpl.getNumScalarDOF();
        }
        
        DOF(const DOF &toCopy) : DOFBase<DataType, PropertyIndex>(toCopy), m_dofImpl() {
            m_dofImpl = toCopy.m_dofImpl;
        }
        
        ~DOF() {
            
        }
        
        constexpr unsigned int getIndex() { return PropertyIndex; }
        
        //The t parameter is optional and is used for spacetime degrees of freedom
        /*inline unsigned int getNumScalarDOF() const {
            return m_dofImpl.getNumScalarDOF();
        }
        
        inline unsigned int getLocalId() const {
            return m_localId;
        }
        
        inline unsigned int getGlobalId() const {
            return m_globalId;
        }
        
        inline void offsetGlobalId(unsigned int offset) {
            m_globalId += offset;
        }*/
        
        //inline std::pair<DataType *, unsigned int> getPtr(const State<DataType> &state) const {
          //  return m_dofImpl.getPtr(state, DOFBase<DataType, PropertyIndex>::m_globalId);
        //}
        
    protected:
        
        Impl<DataType, PropertyIndex> m_dofImpl;
        //unsigned int m_globalId;
        //unsigned int m_localId;
        //unsigned int m_numScalarDOF;
        
    private:
    };
    
}

#endif
