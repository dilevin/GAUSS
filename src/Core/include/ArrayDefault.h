#ifndef _ArrayDefault_H
#define _ArrayDefault_H

#include "CoreDefines.h"
#include "Array.h"

#define DYNAMIC_SIZE_ARRAY 9999

/**
* ArrayDefault - simple main-memory array that is resizable and non-mappable
*/
namespace Core
{
    template<typename TYPE, unsigned int M>
    class ArrayDefault : public Array<TYPE>
	{
        
	public:
        
        explicit ArrayDefault();
        explicit ArrayDefault(const ArrayDefault &toCopy);
        explicit ArrayDefault(const ArrayDefault *toCopy, const unsigned int index);
        
        virtual ~ArrayDefault();
        
		//Accessors
		inline const TYPE & get(const unsigned int index) const;
		inline void set(const unsigned int index, const TYPE &value);
        inline void set(const unsigned int index, unsigned int num, const TYPE *data);
        
		inline TYPE & set(const unsigned int index);
        
        //Special Default Array only methods
        inline TYPE * getPtr();
        inline const TYPE * getConstPtr() const;
        
		//operators
		inline TYPE & operator[](unsigned int index);
        
        //Methods
        ArrayDefault<TYPE,M> * clone() const;
        ArrayDefault<TYPE,M> * map(const unsigned int index = 0) const;
        int resize(size_t newSize);
        
	protected:
        
		TYPE m_memPtr[M]; //ptr to memory
        
	private:
        
		
	};

	template<typename TYPE>
    class ArrayDefault<TYPE, DYNAMIC_SIZE_ARRAY> : public Array<TYPE>
	{

	public:

        explicit ArrayDefault(const size_t numElements=0, const size_t bufferSize  = 0);
        explicit ArrayDefault(const ArrayDefault &toCopy);
        explicit ArrayDefault(const ArrayDefault *toCopy, const unsigned int index);

        virtual ~ArrayDefault();

		//Accessors
		inline const TYPE & get(const unsigned int index) const;
		inline void set(const unsigned int index, const TYPE &value);
        inline void set(const unsigned int index, unsigned int num, const TYPE *data);
		inline TYPE & set(const unsigned int index);

        //Special Default Array only methods
        inline TYPE * getPtr() const;
        inline const TYPE * getConstPtr() const;
        
		//operators
		inline TYPE & operator[](unsigned int index);

        //Methods
        ArrayDefault<TYPE,DYNAMIC_SIZE_ARRAY> * clone() const;
        ArrayDefault<TYPE,DYNAMIC_SIZE_ARRAY> * map(const unsigned int index = 0) const;
        int resize(size_t newSize);

	protected:

		TYPE *m_memPtr; //ptr to memory

	private:

		
	};

}


//Implementation 
namespace Core 
{

    //Dynamic size implementation
	template<typename TYPE>
    ArrayDefault<TYPE,DYNAMIC_SIZE_ARRAY>::ArrayDefault(const size_t numElements, const size_t bufferSize) : Array<TYPE>(numElements, bufferSize)
	{
        this->m_memPtr = new TYPE[this->m_bufferSize];
        this->m_resizable = true;
        this->m_allocated = true;
        this->m_type = GetArrayType<ArrayDefault<TYPE,DYNAMIC_SIZE_ARRAY> >::ArrayType();
        MEMSET(this->m_memPtr, 0, sizeof(TYPE)*this->m_bufferSize);
	}

    template<typename TYPE>
    ArrayDefault<TYPE,DYNAMIC_SIZE_ARRAY>::ArrayDefault(const ArrayDefault &toCopy) : Array<TYPE>(toCopy)
	{

        this->m_bufferSize = toCopy.m_bufferSize;
        this->m_numElements = toCopy.m_numElements;
        this->m_memPtr = toCopy.m_memPtr;
        this->m_resizable = toCopy.m_resizable;
        this->m_mapped = toCopy.m_mapped;

        if(!this->m_mapped) {
            this->m_memPtr = new TYPE[this->m_bufferSize];
            MEMSET(this->m_memPtr, 0, sizeof(TYPE)*this->m_bufferSize);
            MEMCPY(this->m_memPtr, sizeof(TYPE)*this->m_bufferSize, toCopy.m_memPtr, sizeof(TYPE)*this->m_numElements);
        }

		
	}

    template<typename TYPE>
    ArrayDefault<TYPE,DYNAMIC_SIZE_ARRAY>::ArrayDefault(const ArrayDefault *toCopy, const unsigned int index) : Array<TYPE>(toCopy, index)
    {

        //setup pointer here
        this->m_memPtr = &toCopy->m_memPtr[index];
    }

	template<typename TYPE>
    ArrayDefault<TYPE,DYNAMIC_SIZE_ARRAY>::~ArrayDefault()
	{
        if(!this->m_mapped)
            delete[] m_memPtr;
	}

	//Accesors
	template<typename TYPE>
    inline const TYPE & ArrayDefault<TYPE,DYNAMIC_SIZE_ARRAY>::get(const unsigned int index) const
	{
        assert(index < this->m_bufferSize);
		return m_memPtr[index];
	}

	template<typename TYPE>
    inline void ArrayDefault<TYPE,DYNAMIC_SIZE_ARRAY>::set(const unsigned int index, const TYPE &value)
	{
        assert(index < this->m_bufferSize);
		m_memPtr[index] = value;
	}

	template<typename TYPE>
    inline TYPE & ArrayDefault<TYPE,DYNAMIC_SIZE_ARRAY>::set(const unsigned int index)
	{
        assert(index < this->m_bufferSize);
		return m_memPtr[index];
	}

    template<typename TYPE>
    inline void ArrayDefault<TYPE,DYNAMIC_SIZE_ARRAY>::set(const unsigned int index, unsigned int num, const TYPE *data)
    {
        //copy everything quickly
        resize(index+num); //ensure space is available
        std::memcpy((void *)(m_memPtr+index), data, sizeof(TYPE)*num);
    }
	//operators
	template<typename TYPE>
    inline TYPE & ArrayDefault<TYPE,DYNAMIC_SIZE_ARRAY>::operator[](unsigned int index)
	{
        assert(index < this->m_bufferSize);
        return m_memPtr[index];
    }

    template<typename TYPE>
    inline TYPE * ArrayDefault<TYPE,DYNAMIC_SIZE_ARRAY>::getPtr() const
    {
        return m_memPtr;
    }

    template<typename TYPE>
    inline const TYPE * ArrayDefault<TYPE,DYNAMIC_SIZE_ARRAY>::getConstPtr() const
    {
        return m_memPtr;
    }

    
    //Methods
    template<typename TYPE>
    ArrayDefault<TYPE,DYNAMIC_SIZE_ARRAY> * ArrayDefault<TYPE,DYNAMIC_SIZE_ARRAY>::clone() const
    {
        /*ArrayDefault<TYPE> *newArray = new ArrayDefault<TYPE>(this->m_numElements, this->m_bufferSize);

        //copy everything over
        Core::copy(newArray, this);*/

        return NULL;

    }

    template<typename TYPE>
    ArrayDefault<TYPE,DYNAMIC_SIZE_ARRAY> * ArrayDefault<TYPE,DYNAMIC_SIZE_ARRAY>::map(const unsigned int index) const
    {
        return new Core::ArrayDefault<TYPE,DYNAMIC_SIZE_ARRAY>(this, index);
    }


    template<typename TYPE>
    int ArrayDefault<TYPE,DYNAMIC_SIZE_ARRAY>::resize(size_t newSize)
    {
        if(newSize > this->m_bufferSize)
        {
            TYPE *temp = this->m_memPtr;
            this->m_memPtr = new TYPE[newSize];
            this->m_bufferSize = newSize;
            if(this->m_bufferSize > 0)
                MEMCPY(this->m_memPtr, sizeof(TYPE)*this->getSize(), temp, sizeof(TYPE)*this->m_numElements);

            this->m_numElements = newSize;
            delete[] temp;
            return 1;
        }
        //this->m_numElements = newSize;

        return 1;

    }

    //Fixed size implementations
	template<typename TYPE, unsigned int M>
    ArrayDefault<TYPE,M>::ArrayDefault() : Array<TYPE>(M, M)
	{
        this->m_resizable = false;
        this->m_allocated = true;
        this->m_type = GetArrayType<ArrayDefault<TYPE,M> >::ArrayType();
        
        MEMSET(this->m_memPtr, 0, sizeof(TYPE)*M);
	}
    
    template<typename TYPE, unsigned int M>
    ArrayDefault<TYPE,M>::ArrayDefault(const ArrayDefault &toCopy) : Array<TYPE>(toCopy)
	{
        
        this->m_bufferSize = toCopy.m_bufferSize;
        this->m_numElements = toCopy.m_numElements;
        this->m_resizable = false;
        this->m_mapped = false;
        
        MEMCPY(this->m_memPtr, sizeof(TYPE)*M, toCopy.m_memPtr, sizeof(TYPE)*M);
        
	}
    
	template<typename TYPE, unsigned int M>
    ArrayDefault<TYPE,M>::ArrayDefault(const ArrayDefault *toCopy, const unsigned int index) : Array<TYPE>(toCopy, index)
    {
        DEBUGPRINT("Mapping copy constructor not implemented for fixed size arrays\n");
    }
    
	template<typename TYPE, unsigned int M>
    ArrayDefault<TYPE,M>::~ArrayDefault()
	{
	}
    
	//Accesors
	template<typename TYPE, unsigned int M>
    inline const TYPE & ArrayDefault<TYPE,M>::get(const unsigned int index) const
	{
        assert(index < M);
		return m_memPtr[index];
	}
    
	template<typename TYPE,unsigned int M>
    inline void ArrayDefault<TYPE,M>::set(const unsigned int index, const TYPE &value)
	{
        assert(index < M);
		m_memPtr[index] = value;
	}
    
	template<typename TYPE, unsigned int M>
    inline TYPE & ArrayDefault<TYPE,M>::set(const unsigned int index)
	{
        assert(index < M);
		return m_memPtr[index];
	}
    
	//operators
	template<typename TYPE, unsigned int M>
    inline TYPE & ArrayDefault<TYPE,M>::operator[](unsigned int index)
	{
        assert(index < M);
        return m_memPtr[index];
    }
    
	template<typename TYPE, unsigned int M>
    inline TYPE * ArrayDefault<TYPE,M>::getPtr()
    {
        return m_memPtr;
    }
    
	template<typename TYPE, unsigned int M>
    inline const TYPE * ArrayDefault<TYPE,M>::getConstPtr() const
    {
        return m_memPtr;
    }
    
    
    //Methods
	template<typename TYPE, unsigned int M>
    ArrayDefault<TYPE,M> * ArrayDefault<TYPE,M>::clone() const
    {
        DEBUGPRINT("You cannot clone a fixed-size array");
        return NULL;
    }
    
	template<typename TYPE, unsigned int M>
    ArrayDefault<TYPE,M> * ArrayDefault<TYPE,M>::map(const unsigned int /*index*/) const
    {
        DEBUGPRINT("You cannot map a fixed-size array");
        return NULL; //you can't map a stack array
    }
    
    
	template<typename TYPE, unsigned int M>
    int ArrayDefault<TYPE,M>::resize(size_t /*newSize*/)
    {
        
        DEBUGPRINT("You cannot resize a fixed-size array by definition");
        return 0;
        
    }
    
    //Useful typedefs
    typedef ArrayDefault<double, DYNAMIC_SIZE_ARRAY> ArrayDefaultXd; 
}



#endif
