#ifndef _FORCE_H
#define _FORCE_H

namespace Gauss {
    template<typename DataType, typename Impl >
    class Force
    {
    public:
        template<typename ...Params>
        Force(Params ... params) : m_impl(params...) {}
        ~Force() { }
    
        //forces can give you the energy stored, the force itself and the hessian
        template<typename Scalar>
        inline void getEnergy(Scalar &e,  State<DataType> &state) {
            m_impl.energy(e, state);
        }
        
        //forces always act on at least one DOF of the system this function returns which DOF the are acting on.
        auto & getDOF(unsigned int index) {
            return m_impl.getDOF(index);
        }
        
        unsigned int getNumDOF() { return m_impl.getNumDOF(); }
        
        template <typename Vector>
        inline void getForce(Vector &f,  State<DataType> &state) {
            
            m_impl.getForce(f, state);
        }
        
        
        template <typename Matrix>
        inline void getStiffnessMatrix(Matrix &H,  State<DataType> &state) {
            m_impl.getStiffnessMatrix(H, state);
        }
        
        //dummy method which makes forces fit into the standard interface for the assembler
        template<typename Matrix>
        inline void getMassMatrix(Matrix &H,  State<DataType> &state) {
            
        }
        
    protected:
        
        Impl m_impl;
    
    private:
    };
}

#endif

