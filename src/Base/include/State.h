#ifndef STATE_H
#define STATE_H

#include <tuple>
#include "ArrayDefault.h"

namespace Gauss {
    
    /**
     * The state object stores the current state which is made up of all degrees of freedom
     globalId = n for position and velocity
     
     state getVlocitye (n + offset) should be handled in the state using some sort of coordinate term
     
     i.e q queries state(coord 1)
     qDot queries state(coord 2)
     This assumes only 2 variable tyoes (positions and velocities, need to extend for multiple strides)

     */
    template<typename DataType>
    class State
    {
     
    public:
        using StateDataType = DataType;
        
        State(unsigned int offset=0, unsigned int numScalarDOF = 0);
        ~State();
        State(const State &toCopy);
        
        inline unsigned int stateSize(unsigned int i) const {
            return static_cast<unsigned int>((1-2*i)*static_cast<int>(m_offset) +
                                             static_cast<int>(i*m_backingStore.getSize()));
        }
        
        inline unsigned int getNumScalarDOF() { return m_backingStore.getSize(); }
        
        //add scalar DOF and get back index into State
        unsigned int addScalarDOFs(unsigned int numToAdd);
        bool resize(unsigned int newSize);
        
        
        //getters and setters I like.
        inline void setOffset(unsigned int offset) { m_offset = offset; }
        
        //ptr into memory
        std::tuple<DataType *, unsigned int> getStatePtr(unsigned int index = 0) const;
        
        //ptr to a particular category of state variable (either Q or QDot
        template<unsigned int i=0>
        std::tuple<DataType *, unsigned int> getStatePtr(unsigned int offset = 0) const;
        
        //single value
        template<unsigned int i>
        inline DataType & operator[](unsigned int globalId);
        
        inline DataType & operator[] (unsigned int globalId) {
            return m_backingStore[globalId];
        }
        
        //write into and read from backing store
        bool setState(unsigned int index, const DataType * const val, unsigned int numData);
        bool incrementState(unsigned int index, const DataType *val, unsigned int numData);
        
        //property-wise writes and increments
        template<unsigned int i>
        bool setState(const DataType * const val, unsigned int numData);
        
        template<unsigned int i>
        bool incrementState(const DataType *val, unsigned int numData);
        
        
        
    protected:
        
        unsigned int m_offset;
        Core::ArrayDefault<DataType,DYNAMIC_SIZE_ARRAY> m_backingStore;
        
    private:
    };

    template<typename DataType>
    State<DataType>::State(unsigned int offset, unsigned int numScalarDOF) {
        m_offset = offset;
        m_backingStore.resize(numScalarDOF);
    }


    template<typename DataType>
    State<DataType>::~State() {
        
    }

    template<typename DataType>
    State<DataType>::State(const State &toCopy) {
        
        //I'm pretty certain the copy constructor is working
        toCopy.m_backingStore = m_backingStore;
    }

    //add scalar DOF and get back index into State
    template<typename DataType>
    unsigned int  State<DataType>::addScalarDOFs(unsigned int numToAdd) {
        m_backingStore.resize(m_backingStore.getSize() + numToAdd);
        return m_backingStore.getSize();
    }
    
    //resize the state
    template<typename DataType>
    bool  State<DataType>::resize(unsigned int newSize) {
        m_backingStore.resize(newSize);
        return true;
    }
    


    //write into and read from backing store
    template<typename DataType>
    bool State<DataType>::setState(unsigned int index, const DataType * const val, unsigned int numData) {
        //copy data directly into the backing store
        m_backingStore.set(index, numData, val);
    }

    template<typename DataType>
    bool State<DataType>::incrementState(unsigned int index, const DataType *val, unsigned int numData) {
        //increment all values in backing store
        
        assert(index+numData < m_backingStore.getNumElements());
        //FOR LOOP but should be able to parallelize this
        for(unsigned int ii=index; ii < index+numData; ++ii) {
            m_backingStore[ii] += val[ii-index];
        }
    }

    //get ptr into state
    //if
   template<typename DataType>
    std::tuple<DataType *, unsigned int>  State<DataType>::getStatePtr(unsigned int index)  const {
        //get a pointer to the backing store
        
        //make sure you aren't running of the end of the array
        assert(m_backingStore.getSize() > index);
        
        return std::tuple<DataType *, unsigned int>(m_backingStore.getPtr() + index, m_backingStore.getSize()-index);
    }
    
    template<typename DataType>
    template<unsigned int i>
    std::tuple<DataType *, unsigned int>  State<DataType>::getStatePtr(unsigned int offset)  const {
        //get a pointer to the backing store
        
        //make sure you aren't running of the end of the array
        assert(m_backingStore.getSize() > i*m_offset + offset);
        
        return std::tuple<DataType *, unsigned int>(m_backingStore.getPtr() + (m_offset*i+offset),
                                                         stateSize(i) - offset);
    }
    
    
    template<typename DataType>
    template<unsigned int i>
    DataType & State<DataType>::operator[](unsigned int globalId) {
        assert(m_backingStore.getNumElements() > globalId);
        return m_backingStore[globalId + m_offset];
    }
    
    template<typename DataType>
    template<unsigned int i>
    bool State<DataType>::setState(const DataType * const val, unsigned int numData) {
        std::cout<<"State Size:" <<stateSize(i)<<"\n";
        std::cout<<"Backing Size:" <<m_backingStore.getSize()<<"\n";
        assert(numData == stateSize(i));
        
        m_backingStore.set(i*m_offset, numData, val);
        return true;
        
    }
    
    template<typename DataType>
    template<unsigned int i>
    bool State<DataType>::incrementState(const DataType *val, unsigned int numData) {
        assert(numData < stateSize(i));
        //FOR LOOP but should be able to parallelize this
        for(unsigned int ii=i*m_offset; ii < m_offset+stateSize(i); ++ii) {
            m_backingStore[ii] += val[ii-i*m_offset];
        }

    }
    

}

#endif
