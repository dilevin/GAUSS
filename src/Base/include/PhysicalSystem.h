#ifndef _PHYSICALSYSTEM_H
#define _PHYSICALSYSTEM_H

#include "State.h"
#include "MultiVector.h"
#include "Utilities.h"
#include "Assembler.h"

namespace Gauss {
    
    //Needs to be compatible with states of different type
    template<typename DataType, typename Impl>
    class PhysicalSystem
    {
    public:
        
        /*
         *
         * No interface be
         */
        template<typename ...Params>
        PhysicalSystem(Params ...params) : m_systemImpl(params...) { }
        
        PhysicalSystem() : m_systemImpl() { }
        
        ~PhysicalSystem() {
            
        }
        
        template<typename Assembler>
        void getMassMatrix(Assembler &assembler, const State<DataType> &state) const {
            m_systemImpl.getMassMatrix(assembler, state);
        }
        
        template<typename Assembler>
        void getStiffnessMatrix(Assembler &assembler, const State<DataType> &state) {
            m_systemImpl.getStiffnessMatrix(assembler, state);
        }
        
        template <typename Vector>
        inline void getForce(Vector &f, const State<DataType> &state) {
            
            m_systemImpl.getForce(f, state);
        }

        template<typename Assembler>
        void setDOFGlobalIndex(Assembler &assembler) {
            m_systemImpl.setDOFGlobalIndex(assembler);
        }
        
        const auto & getImpl() const { return m_systemImpl; }
        auto & getImpl() { return m_systemImpl; }
        
        auto & getQ() { return m_systemImpl.getQ(); }
        auto & getQ() const { return m_systemImpl.getQ(); }
        
        auto & getQDot() { return m_systemImpl.getQDot(); }
        auto & getQDot() const { return m_systemImpl.getQDot(); }
        
    protected:
        
        Impl m_systemImpl; //system implementation class.
        
    private:
    };
    
    //conveniance function for getting at implementations from pointers
    template<typename DataType, typename Impl>
    inline auto & impl(PhysicalSystem<DataType, Impl> *system) { return system->getImpl(); }
    
}
#endif
