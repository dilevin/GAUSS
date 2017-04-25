#ifndef _ARRAY_H
#define _ARRAY_H

#include "GetArrayType.h"
#include <iostream>
#include <numeric>

/**
* Array - base class for methods that allocate memory on a heap. Array's may be able to dynamically resize themselves, the can be either copied from
* (new allocation) or mapped from (no new allocation just use pointer from other array) an existing array. 
* Copy constructor and the clone method copy the array and contents
* Map produces a mapped array from the given array
* Assignment operator is a wrapper for clone
* All arrays should zero themselves
*/
namespace Core
{
	template<typename TYPE>
	class Array 
	{

	public:

        explicit Array(const size_t numElements = 0, const size_t bufferSize  = 1);
        Array(const Array &toCopy);
        Array(const Array *toCopy, const unsigned int index = 0);

		virtual ~Array();

		//Accesors
        inline size_t getAllocated() const { return m_bufferSize; }
        inline size_t getSize() const { return m_numElements; }
        inline unsigned int getType() const { return m_type; }

		virtual const TYPE & get(const unsigned int index) const;
		virtual void set(const unsigned int index, const TYPE &value);
        virtual void set(const unsigned int index, const unsigned int num, const TYPE *data);
		virtual TYPE & set(const unsigned int index);
        
        inline bool isAllocated() const { return m_allocated; }

		//operators
		virtual TYPE & operator[](unsigned int index);

        //Methods
        virtual Array<TYPE> * clone() const = 0;
        virtual Array<TYPE> * map(const unsigned int index = 0) const = 0;
        virtual int resize(size_t newSize) = 0;

	protected:

		bool m_mapped; //true if array is mapped from another array, false otherwise
		bool m_resizable; //true if array is dynamically resizable, false otherwise
        bool m_allocated;

		size_t m_numElements;
		size_t m_bufferSize;

        unsigned int m_type;

	private:
		
		TYPE failure;
	};

}


//Implementation 
namespace Core 
{
	template<typename TYPE>
    Array<TYPE>::Array(const size_t numElements, const size_t bufferSize)
	{
        //Required by all array types
       // m_type = GetArrayType<Array<TYPE> >::ArrayType();

		//by default arrays are not resizable and not mapped
        m_mapped = false;
        m_resizable = false;

        if(bufferSize == 0) {
            // Set bufferSize=numElements by default
            m_bufferSize = numElements;
        } else if(bufferSize < numElements) {
			std::cout<<"bufferSize < numElements: Setting bufferSize = numElements \n";
			m_bufferSize = std::max<size_t>(numElements,1); //Robustly handle initializing empty arrays
		} else {
			m_bufferSize = bufferSize;
		}

		m_numElements = numElements;
        m_allocated = false; //nothing allocated yet
	}

    template<typename TYPE>
    Array<TYPE>::Array(const Array *toCopy, const unsigned int index)
    {
        m_mapped = true;
        m_resizable = false; //Mapped arrays are fixed
        m_numElements = toCopy->m_numElements - index;
        m_bufferSize = toCopy->m_bufferSize - index;
        m_type = toCopy->m_type;
    }

	template<typename TYPE>
    Array<TYPE>::Array(const Array &/*toCopy*/)
	{
		std::cout<<"Array base class copy constructor \n";
	}

	template<typename TYPE>
	Array<TYPE>::~Array()
	{
		//Default destructor does nothing
	}

	//Accesors
	template<typename TYPE>
    const TYPE & Array<TYPE>::get(const unsigned int /*index*/) const
	{
		std::cout<<"Array base class get method \n";
		return failure; //This does nothing, just gets around the compiler complaining that I need to return something
	}

	template<typename TYPE>
    void Array<TYPE>::set(const unsigned int /*index*/, const TYPE &/*value*/)
	{
		std::cout<<"Array base class set method \n";
	}

	template<typename TYPE>
    TYPE & Array<TYPE>::set(const unsigned int /*index*/)
	{

		std::cout<<"Array base class set method \n";
		return failure;
	}

    template<typename TYPE>
    void Array<TYPE>::set(const unsigned int index, unsigned int num, const TYPE *data)
    {
        std::cout<<"Array base class batch set method \n";
    }
    
	//operators
	template<typename TYPE>
    TYPE & Array<TYPE>::operator[](unsigned int /*index*/)
	{
		std::cout<<"Array base class [] operator \n";
		return failure;
	}

}

#endif
