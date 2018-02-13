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
        
        constexpr static unsigned int numTypes() { return sizeof...(Types); }

        template <int... M>
        using IS = typename std::integer_sequence<int,M...>;
        constexpr static auto _is() { return std::make_integer_sequence<int,numTypes()>(); }

        unsigned int getNumCategories() const { return numTypes();}
        
        template<unsigned int systemType>
        inline auto&& get() {
            static_assert(systemType < numTypes(),"");
            
            return std::get<systemType>(m_vectorTuple);
        }
        template<unsigned int systemType>
        inline auto get(unsigned int index) {
            auto&& myvec = get<systemType>();
            assert(index < myvec.size());
            
            return myvec[index];
        }
    template <int... M, typename Func>
    void forEachVector(IS<M...>, Func&& f) {
#ifdef __cpp_fold_expressions
        (f(get<M>()),...); 
#else
        (void)std::initializer_list<int>{ (f(get<M>()),0)... };
#endif
    }
    template <typename Func>
    void forEachVector(Func&& f) { forEachVector(_is(),std::forward<Func>(f)); }

    template <int... M, typename Func>
    void forEachVectorWithClassIndex(IS<M...>, Func&& f) {
#ifdef __cpp_fold_expressions
        (f(M,get<M>()),...); 
#else
        (void)std::initializer_list<int>{(f(M,get<M>()),0)...};
#endif
    }
    template <typename Func>
    void forEachVectorWithClassIndex(Func&& f) { forEachVectorWithClassIndex(_is(),std::forward<Func>(f)); }

    template <typename Func>
    void forEach(Func&& f) { forEachVector([&](auto&& vec) {
            for(auto&& v: vec) {
                f(v);
            }
            }); 
    }
    template <typename Func>
    void forEachIndex(Func&& f) { forEachVectorWithClassIndex([&](int ci, auto&& vec) {
            for(unsigned int i = 0; i < vec.size(); ++i) {
                f(ci,i,vec[i]);
            }
            }); 
    }



        
    protected:
           std::tuple<std::vector<Types>...> m_vectorTuple;
        
    private:
    };
    
    
    //utility methods

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
    

    template<typename ...T, typename Func>
    inline void forEach(MultiVector<T...> &mv, Func func) {
        mv.forEach(func);
    }

    template<typename ...T, typename Func>
    inline void forEachIndex(MultiVector<T...> &mv, Func func) {
        mv.forEachIndex(func);
    }
    
    
    template<typename ...T, typename Func>
    inline void eachVector(MultiVector<T...> &mv, Func func) {
        mv.forEachVector(func);
    }
    
    //member functions that need utility code above
    template<typename ...Types>
    inline void MultiVector<Types...>::clear() {
        eachVector(*this, MEMBER_FUNC(clear, )); //empty space for a noparam
    }
    
}

#endif /* MultiVector_h */
