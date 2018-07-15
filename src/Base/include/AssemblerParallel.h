//
//  AssemblerParrallel.h
//  Gauss
//
//  Created by David Levin on 6/6/17.
//
//

#ifndef AssemblerParallel_h
#define AssemblerParallel_h

#ifdef GAUSS_OPENMP

#include <omp.h>
#include <Assembler.h>
#include <CoreDefines.h>
#include <Utilities.h>
#include <UtilitiesOMP.h>
#include <World.h>

namespace Gauss {
    
    
    //A parallel assembler just uses one serial assembler per available thread
    template<typename SerialAssembler>
    class AssemblerParallelImpl : public AssemblerBase {
    public:
        
        using MatrixType = typename SerialAssembler::MatrixType;
        
        AssemblerParallelImpl() {
            
            //Number of available theads
            std::cout<<"Number of Available Threads: "<<omp_thread_count()<<"\n";
            
            
            m_serialAssemblers.resize(omp_thread_count());
            
        }
        
        inline void init(unsigned int m, unsigned int n=1, unsigned int rowOffset = 0, unsigned int colOffset = 0) {
         
            //do everything in parallel
            m_assembled.resize(m,n); // TODO need conditional to deal with case where MatrixType is a vector
            m_assembled.setZero();
            
            for(unsigned int ii=0; ii < m_serialAssemblers.size(); ++ii) {
                m_serialAssemblers[ii].init(m,n, rowOffset, colOffset);
            }
        }
        
        
        inline void finalize() {
            
            //build giant triplets list and set it up
            //I think I want to assemble seperately then add (split up my setTriplets time)
            #pragma omp parallel
            {
                #pragma omp for
                for(unsigned int ii=0; ii < m_serialAssemblers.size(); ++ii) {
                    m_serialAssemblers[ii].finalize();
                }
            }
    
            for(unsigned int ii=0; ii<m_serialAssemblers.size(); ++ii) {
                m_assembled += (*m_serialAssemblers[ii]);
            }
            
        }
        
        //convenient overloads
        inline auto & getMatrix() {
            return m_assembled;
            
        }
        
        template<typename I, typename J, typename Input>
        inline void assemble(I &i, J &j, Input &toAssembler) {
            m_serialAssemblers[0].getImpl().assemble(i,j, toAssembler); //default single threaded behavior
            //exit(1);
        }
        
        template<typename I, typename Input>
        inline void assemble(I &i, Input &toAssembler) {
            m_serialAssemblers[0].getImpl().assemble(i, toAssembler); //default single threaded behavior
            //exit(1);
        }
        
        //next step, this needs to change to take in a list of i's, j's and sizes
        //take in std::vectors for indices and size
        
        inline void setOffset(unsigned int rowOffset, unsigned int colOffset = 0) {
            AssemblerBase::setOffset(rowOffset, colOffset);
            
            #pragma omp parallel
            {
                #pragma omp for
                for(unsigned int ii=0; ii < m_serialAssemblers.size(); ++ii) {
                    m_serialAssemblers[ii].setOffset(rowOffset, colOffset);
                }
                
            }
        }
    
        inline SerialAssembler & operator[](unsigned int threadId) {
            return m_serialAssemblers[threadId];
        }
        
        SerialAssembler & getAssembler(unsigned int threadId) {
            return m_serialAssemblers[threadId];
        }
        
        
        //For MVP assemblers
        inline void setX(MatrixType &x) {
        #pragma omp parallel
            {
            #pragma omp for
                for(unsigned int ii=0; ii < m_serialAssemblers.size(); ++ii) {
                    (m_serialAssemblers[ii].getImpl().setX(x));
                }
            }
            
            return *this;
        }
        
        //handle operators
        template<typename Params>
        inline AssemblerParallelImpl & operator*=(Params &x) {
            #pragma omp parallel
            {
                #pragma omp for
                for(unsigned int ii=0; ii < m_serialAssemblers.size(); ++ii) {
                    (m_serialAssemblers[ii].getImpl())*=x;
                }
           }
            
            return *this;
        }
        
        
    protected:
        
        std::vector<SerialAssembler> m_serialAssemblers;
        
        //At some point I should figure out this type based on the constituent asemblers but for now just assume Eigen
        typename SerialAssembler::ImplType::MatrixType m_assembled;
        
    private:
    };
    
    template<typename DataType, typename SerialAssembler>
    using AssemblerParallel = Assembler<DataType, AssemblerParallelImpl<SerialAssembler> >;
 
    template<typename DataType, typename SerialAssembler>
    struct IsParallel<Assembler<DataType, AssemblerParallelImpl<SerialAssembler> > > {
    public:
        constexpr static bool value = true;
    };
}
#else
    template<typename DataType, typename SerialAssembler>
    using AssemblerParallel = SerialAssembler;
#endif //OPENMP is Available

#endif /* AssemblerParrallel_h */
    
    
