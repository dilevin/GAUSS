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
        
        using ImplType = Impl;
        /*
         *
         * No interface be
         */
        template<typename ...Params>
        PhysicalSystem(Params ...params) : m_systemImpl(params...) { }
        
        PhysicalSystem() : m_systemImpl() { }
        
        ~PhysicalSystem() {
            
        }
        
        inline DataType getEnergy(const State<DataType> &state) { return m_systemImpl.getEnergy(state); }
        inline DataType getStrainEnergy(const State<DataType> &state) { return m_systemImpl.getStrainEnergy(state); }
        inline decltype(auto) getStrainEnergyPerElement(const State<DataType> &state) { return m_systemImpl.getStrainEnergyPerElement(state); }
        
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

        template <typename Vector>
        inline void getInternalForce(Vector &f, const State<DataType> &state) {
            m_systemImpl.getInternalForce(f, state);
        }

        template<typename Assembler>
        void setDOFGlobalIndex(Assembler &assembler) {
            m_systemImpl.setDOFGlobalIndex(assembler);
        }
        
        template<typename ...Params>
        inline decltype(auto) getPosition(const State<DataType> &state, Params &...params) const {
            return m_systemImpl.getPosition(state, params...);
        }
        
        template<typename Vector, typename ...Params>
        inline auto getDPDQ(Vector &x, Params &... params) {
            return m_systemImpl.getDPDQ(x, params...);
        }
        
        template<typename  Vector, typename ...Params>
        inline auto getVelocity(Vector &x, Params &...params) {
            return m_systemImpl.getVelocity(x, params...);
        }
        
        template<typename  Vector, typename ...Params>
        inline auto getDVDQ(Vector &x, Params &...params) {
            return m_systemImpl.getDVDQ(x, params...);
        }
        
        inline const auto & getImpl() const { return m_systemImpl; }
        inline auto & getImpl() { return m_systemImpl; }
        
        //need more of these so i can get the Q's that support a spatial point
        template<typename ...Params>
        inline decltype(auto) getQ(Params &...params) { return m_systemImpl.getQ(params...); }
        
        template<typename ...Params>
        inline decltype(auto) getQ(Params &...params) const { return m_systemImpl.getQ(params...); }
        
        template<typename ...Params>
        inline decltype(auto) getQDot(Params &...params) { return m_systemImpl.getQDot(params...); }
        
        template<typename ...Params>
        inline decltype(auto) getQDot(Params &...params) const { return m_systemImpl.getQDot(params...); }
        
        
        //get geometry (had to add this for collision detection code)
        inline auto getGeometry() { return m_systemImpl.getGeometry(); }
        
    protected:
        
        Impl m_systemImpl; //system implementation class.
        
    private:
    };
    
    //conveniance function for getting at implementations from pointers
    template<typename DataType, typename Impl>
    inline auto & impl(PhysicalSystem<DataType, Impl> *system) { return system->getImpl(); }
    
}
#endif
