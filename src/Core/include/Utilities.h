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
