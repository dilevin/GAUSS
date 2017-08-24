//
//  TimeStepperLinearModes.h
//  Gauss
//
//  Created by David Levin on 8/23/17.
//
//

#ifndef TimeStepperLinearModes_h
#define TimeStepperLinearModes_h

#include <cmath>

//This time stepper applies osscilatory displacements from a set of linear modes
//Just used to display modes
namespace Gauss {
    
    //Given Initial state, step forward in time using linearly implicit Euler Integrator
    template<typename DataType>
    class TimeStepperImplLinearModes
    {
    public:
        
        TimeStepperImplLinearModes(Eigen::MatrixXx<DataType> &modes, Eigen::VectorXx<DataType> &frequencies, unsigned int modeIndex = 6) {
         
            m_modes = modes;
            m_frequencies = frequencies;
            m_modeIndex = modeIndex; //first non-rigid mode
            
        }
        
        TimeStepperImplLinearModes(const TimeStepperImplLinearModes &toCopy) {
            
        }
        
        ~TimeStepperImplLinearModes() { }
        
        //Methods
        //init() //initial conditions will be set at the begining
        template<typename World>
        void step(World &world, double dt, double t);
        
    protected:
        
        Eigen::MatrixXx<DataType> m_modes;
        Eigen::VectorXx<DataType> m_frequencies;
        unsigned int m_modeIndex; //selected mode for display
        
    private:
    };

    template<typename DataType>
    template<typename World>
    void TimeStepperImplLinearModes<DataType>::step(World &world, double dt, double t) {
    
        //make osscilatory animation
        Eigen::Map<Eigen::VectorXd> q = mapStateEigen<0>(world);
        Eigen::Map<Eigen::VectorXd> qDot = mapStateEigen<1>(world);
        
        q = 100.0*m_modes.col(m_modeIndex)*std::cos(t);
        qDot = -m_modes.col(m_modeIndex)*std::sin(t);
    
    }

    
    template<typename DataType>
    using TimeStepperLinearModes = TimeStepper<DataType, TimeStepperImplLinearModes<DataType> >;

}




#endif /* TimeStepperLinearModes_h */
