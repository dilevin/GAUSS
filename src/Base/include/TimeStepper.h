#ifndef _TIMESTEPPER_H
#define _TIMESTEPPER_H

#include <World.h>
#include <Assembler.h>
//TODO Solver Interface
//TODO External Forces
//TODO Constraints
namespace Gauss {
    
    
    //Given Initial state, step forward in time using linearly implicit Euler Integrator
    template<typename DataType, typename Impl>
    class TimeStepper
    {
    public:
        
        TimeStepper(DataType dt) : m_impl() {  m_dt = dt; m_t = 0;}
        
        TimeStepper(const TimeStepper &toCopy) {
            m_dt = toCopy.m_dt;
            m_t = toCopy.m_t;
            m_impl = toCopy.m_impl;
        }
        
        ~TimeStepper() { }
        
        //Methods
        template<typename World>
        void step(World &world) { m_impl.step(world, m_dt); m_t += m_dt; }
        inline DataType getTime() const { return m_t; }
        inline void setDt(DataType dt) { m_dt = dt; }
        inline auto & getLagrangeMultipliers() { return m_impl.getLagrangeMultipliers(); }
        
    protected:
        
        Impl m_impl;
        DataType m_t;
        DataType m_dt;
        
    private:
    };
}

#endif
