//
//  MultiVector.h
//  Gauss
//
//  Created by David Levin on 1/29/17.
//
//

#ifndef MultiVector_h
#define MultiVector_h

#include <iostream>
#include <tuple>
#include <vector>

//Generic Lambda for calling member functions of objects stored in the multivector
#define MEMBER_FUNC(funcName, param) [](auto&& obj) -> decltype(auto) \
                        { return std::forward<decltype(obj)>(obj).funcName(param); } \


namespace Gauss {
    
    template<typename ...Types>
    class MultiVector
    {
    public:
     
        using TupleType = std::tuple<std::vector<Types>...>;
        
        //add
        template<typename T>
        inline void add(T toAdd) { std::get<std::vector<T>>(m_vectorTuple).push_back(toAdd); }
        
        //remove
        inline void clear();
        
        //access
        std::tuple<std::vector<Types>...> & getStorage() { return m_vectorTuple; }
        
        unsigned int getNumCategories() const { return std::tuple_size<std::tuple<std::vector<Types>...> >::value;}
        
        constexpr unsigned int numTypes() { return std::tuple_size<std::tuple<std::vector<Types>...> >::value; }
        
        template<unsigned int systemType>
        inline auto get(unsigned int index) {
            assert(index < numTypes());
            assert(index < std::get<systemType>(m_vectorTuple).size());
            
            return std::get<systemType>(m_vectorTuple)[index];
        }
        
    protected:
           std::tuple<std::vector<Types>...> m_vectorTuple;
        
    private:
    };
    
    
    //utility methods
    template<typename T, typename Func,unsigned int index>
    class forEachClass
    {
    public:
        inline forEachClass(T &tuple, Func &f) {
        
            forEachClass<T, Func, index-1>(tuple, f);
        
            for(auto &itr : std::get<index>(tuple))
            {
                f(itr);
            }
        
        }
    };
    
    //terminate recursion
    template<typename T, typename Func>
    class forEachClass<T,Func,0>
    {
    public:
        inline forEachClass(T &tuple, Func &f) {
            
            for(auto &itr : std::get<0>(tuple))
            {
                f(itr);
            }
            
        }
    };

    //these loops give me the tuple and vector index of the paricular item
    template<typename T, typename Func,unsigned int index>
    class forEachIndexClass
    {
    public:
        inline forEachIndexClass(T &tuple, Func &f) {
            
            forEachIndexClass<T, Func, index-1>(tuple, f);
            
            for(unsigned int ii=0; ii<std::get<index>(tuple).size(); ++ii)
            {
                f(index, ii, std::get<index>(tuple)[ii]);
            }
            
        }
    };
    
    //terminate recursion
    template<typename T, typename Func>
    class forEachIndexClass<T,Func,0>
    {
    public:
        inline forEachIndexClass(const T &tuple, Func &f) {
            
            for(unsigned int ii=0; ii<std::get<0>(tuple).size(); ++ii)
            {
                f(0, ii, std::get<0>(tuple)[ii]);
            }
        }
    };

    //convenience functions
    template<typename ...T, typename Func>
    inline void forEach(MultiVector<T...> &mv, Func func) {
        forEachClass<typename MultiVector<T...>::TupleType,Func,std::tuple_size<typename MultiVector<T...>::TupleType>::value-1>(mv.getStorage(),func);
    }
    
    template<typename ...T, typename Func>
    inline void forEachIndex(MultiVector<T...> &mv, Func func) {
        forEachIndexClass<typename MultiVector<T...>::TupleType,Func,std::tuple_size<typename MultiVector<T...>::TupleType>::value-1>(mv.getStorage(),func);
    }
    
    //recursively touch each vector in the multivector (could probably build all these functions out of
    //lambdas but in a rush now)
    template<typename T, typename Func,unsigned int index>
    class eachVectorClass
    {
    public:
        inline eachVectorClass(T &tuple, Func &f) {
            
            eachVectorClass<T, Func, index-1>(tuple, f);
            f(std::get<index>(tuple));
        }
    };
    
    template<typename T, typename Func>
    class eachVectorClass<T,Func,0>
    {
    public:
        inline eachVectorClass(T &tuple, Func &f) {
            
                f(std::get<0>(tuple));
        }
    };
    
    //convenience functions
    template<typename ...T, typename Func>
    inline void eachVector(MultiVector<T...> &mv, Func func) {
        eachVectorClass<typename MultiVector<T...>::TupleType,Func,std::tuple_size<typename MultiVector<T...>::TupleType>::value-1>(mv.getStorage(),func);
    }
    
    //member functions that need utility code above
    template<typename ...Types>
    inline void MultiVector<Types...>::clear() {
        eachVector(*this, MEMBER_FUNC(clear, )); //empty space for a noparam
    }
    
}

#endif /* MultiVector_h */
