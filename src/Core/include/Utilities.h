//
//  Utilities.h
//  Gauss
//
//  Created by David Levin on 1/31/17.
//
// 

#ifndef Utilities_h
#define Utilities_h

#ifdef GAUSS_OPENMP
#include <omp.h>
#endif

#include "State.h"

//Eigen Stuff
#include <Eigen/Dense>
#include <Eigen/Sparse>

#define STRINGIFY(s) #s

#define DataDir(s) STRINGIFY(s)

#define Vert(f, c) V(F(f), c) //Handle raw vertex, face pointers

#define Vec3f(x,y,z) std::array<float,3>({{x,y,z}}).data()

#define Vec3d(x,y,z) std::array<double,3>({{x,y,z}}).data()

//Random utilities and classes that I might need
namespace Gauss {
    
    //Assignment operators to handle interoperability of differenty types
    
    template<typename Src, typename Dst>
    class increment {
    public:
        inline increment(Src &src, const Dst &dst) { std::cout<<"Default increment operator \n"; }
    protected:
    private:
    };
    
    template<typename Src, typename Dst>
    class update {
    public:
        inline update(Src &src, const Dst &dst) { std::cout<<"Default increment operator \n"; }
    protected:
    private:
    };
    
    //increment specializations
    
    //state and Eigen::Vectors
    
    //static if using templates taken from
    //http://baptiste-wicht.com/posts/2015/07/simulate-static_if-with-c11c14.html
    struct identity {
        template<typename T>
        inline T operator()(T&& x) const {
            return std::forward<T>(x);
        }
    };
    
    template<bool Cond>
    struct statement {
        template<typename F>
        inline void then(const F& f){
            f(identity());
        }
        
        template<typename F>
        inline void else_(const F&){}
    };
    
    template<>
    struct statement<false> {
        template<typename F>
        inline void then(const F&){}
        
        template<typename F>
        inline void else_(const F& f){
            f(identity());
        }
    };
    
    template<bool Cond, typename F>
    inline statement<Cond> static_if(F const& f){
        statement<Cond> if_;
        if_.then(f);
        return if_;
    }
    
    //deal with arrays
    template<typename ...T>
    class ArrayUint {
        
    public:
        inline ArrayUint(T... vals) : m_val{vals...} {
            
        }
        
        inline unsigned int operator[](unsigned int index) { assert(index < m_N); return m_val[index]; }
        constexpr unsigned int size() { return m_N; }
    protected:
        unsigned int m_N = sizeof...(T);
        unsigned int m_val[sizeof...(T)];
    };
    
    std::string dataDir();
    
    std::string timeStampString(std::string toStamp);
    
    //Make handling parallel stuff easier
    template<bool IsParallel>
    struct forLoop {
        template <typename T, typename Assembler, typename Func>
        inline forLoop(T &iterateOver, Assembler &assembler, Func &&f) {
            
            //iterate
            for(auto &itr : iterateOver)
            {
                f(assembler, itr);
            }
        }
    };
    
    //is parrallel checker
    template<typename Obj>
    struct IsParallel {
    public:
        constexpr static bool value = false;
    };
   
    //Direct access into a tuple to run a designated function
	#if defined(_WIN32) || defined(_WIN64) || defined (WIN32)
        //Slow version that doesn't break Visual Studio Compiler
        template<unsigned int CheckIndex>
        class ApplyTuple {
        public:
            template<typename Tuple, typename Func, typename ...Params>
            inline static decltype(auto) apply(Tuple &tuple, unsigned int index, Func &func, Params ...params) {
                if(index == CheckIndex) {
                    return func(std::get<CheckIndex>(tuple), params...);
                }
                
                return ApplyTuple<CheckIndex-1>::apply(tuple,index, func, params...);
            }
            
        };
        
        template<>
        class ApplyTuple<0> {
        public:
            template<typename Tuple, typename Func, typename ...Params>
            inline static decltype(auto) apply(Tuple &tuple, unsigned int index, Func &func, Params ...params) {
                if(index == 0) {
                    return func(std::get<0>(tuple), params...);
                }
                
                std::cout<<"Error in ApplyTuple, no index found \n";
                exit(0);
                
            }
            
        };
        
        template<typename Tuple, typename Func, typename ...Params>
        inline decltype(auto) apply(Tuple &tuple, unsigned int index, Func &func, Params ...params){
            
            return ApplyTuple<std::tuple_size<Tuple>::value -1>::apply(tuple, index, func, params...);
        }
    #else
        //O(1) version for gcc and clang
        template <typename T>
        struct FunctionProperties  {
            using ReturnType = void;
        };
        
        template <typename CT, typename RT, typename... Args>
        struct FunctionProperties<RT(CT::*)(Args...) const > {
            using ReturnType = RT;
        };
        
        template <typename CT, typename RT, typename... Args>
        struct FunctionProperties<RT(CT::*)(Args...)> {
            using ReturnType = RT;
        };
        
        // for function pointers
        template <typename RT, typename... Args>
        struct FunctionProperties<RT (*)(Args...)>  {
            using ReturnType = RT;
        };
        
        //virutal function-like behaviour for tuples
        template<int N, class Tuple, class FunctionWrapper, typename ...Params>
        inline auto apply_one(Tuple & p, FunctionWrapper &func, Params ...params)
        {
            return func(std::get<N>(p), params...);
        }
        
        //define function table
        template<typename A, typename B, typename C,typename D, typename ...E>
        class FunctionTable { };
        
        template<std::size_t... Is, typename Tuple, typename FunctionWrapper, typename ReturnType, typename ...Params>
        class FunctionTable<std::index_sequence<Is...>, Tuple, FunctionWrapper, ReturnType, Params...> {
        public:
            static ReturnType (*lookup_table[std::tuple_size<Tuple>::value])(Tuple&, FunctionWrapper &, Params ...);
        };
        
        template<std::size_t... Is, typename Tuple, typename FunctionWrapper, typename ReturnType, typename ...Params>
        ReturnType (*FunctionTable<std::index_sequence<Is...>, Tuple, FunctionWrapper, ReturnType, Params...>::lookup_table[std::tuple_size<Tuple>::value])(Tuple&, FunctionWrapper &, Params ... )  = {  &apply_one<Is, Tuple, FunctionWrapper, Params...>... };
        
        
        template<typename Tuple, typename Func, typename ...Params>
        inline decltype(auto) apply(Tuple &tuple, unsigned int index, Func &func, Params ...params) {
            using ReturnType = decltype(func(std::get<0>(tuple), params...));
            return FunctionTable<std::make_index_sequence<std::tuple_size<typename std::remove_reference<decltype(tuple)>::type>::value>, typename std::remove_reference<decltype(tuple)>::type, Func, ReturnType, Params... >::lookup_table[index](tuple, func, params...);
        }
    #endif
    
   

    
#ifdef GAUSS_OPENMP
    //Parallel version
    template<>
    struct forLoop<true> {
        template <typename Func, typename Assembler, typename T>
        inline forLoop(T &iterateOver, Assembler &assembler, Func &&f) {
            
            #pragma omp parallel
            {
                #pragma omp for
                //iterate
                for(unsigned int ii=0; ii < iterateOver.size(); ++ii)
                {
                    f(assembler.getImpl()[omp_get_thread_num()], iterateOver[ii]);
                }
            }
        }
        
    };
#endif
    
    
}
#endif /* Utilities_h */
