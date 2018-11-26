#ifndef ASSEMBLER_H
#define ASSEMBLER_H

#include <CoreDefines.h>
#include <Utilities.h>
#include <World.h>

//Needed for more complex operations
#define ASSEMBLEMATINIT(assembler, mSize, nSize) \
assembler.init(mSize, nSize);

#define ASSEMBLEVECINIT(assembler, mSize) \
assembler.init(mSize);

#define ASSEMBLEEND(assembler) \
assembler.finalize();

#define ASSEMBLELIST(assembler, list, funcName) \
forEach(list, [&world, &assembler](auto a) { \
    a->funcName(assembler, world.getState()); \
});

#define ASSEMBLELISTCONSTRAINT(assembler, list, funcName) \
forEach(list, [&world, &assembler](auto a) { \
a->funcName(assembler, world, world.getState()); \
});

#define ASSEMBLELISTOFFSET(assembler, list, funcName, rowOffset, colOffset) \
assembler.setOffset(rowOffset, colOffset);\
forEach(list, [&world, &assembler](auto a) { \
a->funcName(assembler, world, world.getState()); \
});

#define ASSEMBLELISTOFFSETTRANSPOSE(assembler, list, funcName, rowOffset, colOffset) \
assembler.setOffset(rowOffset, colOffset);\
forEach(list, [&world, &assembler](auto a) { \
a->template funcName<decltype(world), decltype(assembler), 1>(assembler, world, world.getState()); \
});

//ASSEMBLEBLAH(list)
//ENDASSEMBLE()
//a = b
/*#define EQUALSMAT(dofsI, dofsJ, a,b) \
static_if<isAssembler<Assembler>::val()>([&](auto f){ \
    f(a).set(dofsI, dofsJ, b);\
}).else_([&](auto f) {\
    f(a) = b;\
});*/

// EQUALS MAT should take the connectivity as a variable
namespace Gauss {
    
    //utility functions for converting actual objects to pointers (smooth out some interface issues)
    template<typename T>
    inline T* ptr(T* obj)  { return obj; }
    
    template<typename T>
    inline T* ptr(T &obj) { return &obj; }
    
    class AssemblerBase {
    public:
        AssemblerBase() { m_rowOffset = 0; m_colOffset = 0; }
        
        inline void setOffset(unsigned int rowOffset, unsigned int colOffset) {
            m_rowOffset = rowOffset;
            m_colOffset = colOffset;
        }
        
    protected:
        int m_rowOffset;
        int m_colOffset;
    private:
        
    };

    //could looping code live here ?
    //ITERATE(list)
    //either a parallel for loop or not
    
    template<typename DataType, typename Impl>
    class Assembler : public AssemblerBase {
    public:
        
        using ImplType = Impl;
        using MatrixType = typename Impl::MatrixType;
        
        Assembler() : m_impl() {  }
        ~Assembler() { }
        
        inline void init(unsigned int m, unsigned int n=1, int rowOffset = 0, int colOffset = 0) {
            m_impl.setOffset(0,0);
            m_impl.init(m,n);
        }
        
        
        inline void finalize() { m_impl.finalize(); }
        
        //convenient overloads
        inline auto & operator*() { return m_impl.getMatrix(); }
        inline const auto & operator*() const { return m_impl.getMatrix(); }
        
        inline auto operator-() { return -m_impl.getMatrix(); }
        inline const auto operator-() const { return -m_impl.getMatrix(); }
        
        
        //next step, this needs to change to take in a list of i's, j's and sizes
        //take in std::vectors for indices and size
        
        //set Transpoe would just swap i and j 
        template<typename I, typename J, typename Input>
        inline void set(I i, J j, Input &toAssemble) {
            m_impl.assemble(i, j, toAssemble);
        }
        
        template<typename I, typename Input>
        inline void set(I i, Input &toAssemble) {
            m_impl.assemble(i, toAssemble);
        }
        
        inline void setOffset(int rowOffset, int colOffset = 0) {
            m_impl.setOffset(rowOffset, colOffset);
        }
        
        Impl & getImpl() { return m_impl; }
        const Impl & getImpl() const { return m_impl; }
    
    private:
    
        Impl m_impl; //specific implementation
        int m_rowOffset;
        int m_colOffset;
    };
}


#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace Gauss {
    class AssemblerImplEigenSparseMatrix : public AssemblerBase {
        typedef double Precision;
        
    public:
        
        using MatrixType = Eigen::SparseMatrix<Precision, Eigen::RowMajor>;
        
        using AssemblerBase::m_rowOffset;
        using AssemblerBase::m_colOffset;
        
        template<typename I, typename J, typename ...Input>
        struct assembleStruct {
            
            inline assembleStruct(AssemblerImplEigenSparseMatrix *parent, I &i, J &j, const double &toAssembler) {
                
                //add double to system as a diagonal matrix
                for(unsigned int idof=0;idof < i.size(); ++idof) {
                    for(unsigned int jdof=0;jdof < j.size(); ++jdof) {
                        for(unsigned int ii=0; ii<ptr(i[idof])->getNumScalarDOF(); ++ii) {
                            parent->m_tripletList.push_back(Eigen::Triplet<Precision>(parent->m_rowOffset+ptr(i[idof])->getGlobalId()+ii,
                                                                                      parent->m_colOffset+ptr(j[jdof])->getGlobalId()+ii,
                                                                                      toAssembler));
                        }
                    }
                }
            }
            
            template<typename Derived>
            inline assembleStruct(AssemblerImplEigenSparseMatrix *parent, I &i, J &j, const Eigen::MatrixBase<Derived> &toAssembler) {
                
                unsigned int ii = 0;
                unsigned int jj = 0;
                unsigned int idof = 0;
                unsigned int jdof = 0;
                unsigned int ipos, jpos;
                
                ipos = 0;
                for(idof=0, ii=0; idof<i.size(); ++idof) {
                    for(ii=0; ii<ptr(i[idof])->getNumScalarDOF(); ++ii) {
                        jpos = 0;
                        for(jdof = 0, jj=0; jdof<j.size(); ++jdof) {
                            for(jj = 0; jj<ptr(j[jdof])->getNumScalarDOF(); ++jj) {
                                parent->m_tripletList.push_back(Eigen::Triplet<Precision>(parent->m_rowOffset+ptr(i[idof])->getGlobalId()+ii,
                                                                                          parent->m_colOffset+ptr(j[jdof])->getGlobalId()+jj,
                                                                                          toAssembler(ipos+ii, jpos+jj)));
                            }
                            
                            jpos += jj;
                        }
                    }
                    
                    ipos += ii;
                }
            }
            
        };
        
        /*template<typename I, typename J>
        struct assembleStruct<I, J, const double> {
            inline assembleStruct(AssemblerImplEigenSparseMatrix *parent, I &i, J &j, const double &toAssembler) {
                
                //add double to system as a diagonal matrix
                for(unsigned int idof=0;idof < i.size(); ++idof) {
                    for(unsigned int jdof=0;jdof < j.size(); ++jdof) {
                            for(unsigned int ii=0; ii<ptr(i[idof])->getNumScalarDOF(); ++ii) {
                                parent->m_tripletList.push_back(Eigen::Triplet<Precision>(parent->m_rowOffset+ptr(i[idof])->getGlobalId()+ii,
                                                                                          parent->m_colOffset+ptr(j[jdof])->getGlobalId()+ii,
                                                                                          toAssembler));
                            }
                    }
                }
            }
        };*/
        
        /*template<typename I, typename J, int ROWS, int COLS>
        struct assembleStruct<I, J, Eigen::Matrix<double, ROWS, COLS> > {
            inline assembleStruct(AssemblerImplEigenSparseMatrix *parent, I &i, J &j, const Eigen::Matrix<double, ROWS, COLS> &toAssembler) {
                
                unsigned int ii = 0;
                unsigned int jj = 0;
                unsigned int idof = 0;
                unsigned int jdof = 0;
                unsigned int ipos, jpos;
            
                ipos = 0;
                for(idof=0, ii=0; idof<i.size(); ++idof) {
                    for(ii=0; ii<ptr(i[idof])->getNumScalarDOF(); ++ii) {
                        jpos = 0;
                        for(jdof = 0, jj=0; jdof<j.size(); ++jdof) {
                            for(jj = 0; jj<ptr(j[jdof])->getNumScalarDOF(); ++jj) {
                                parent->m_tripletList.push_back(Eigen::Triplet<Precision>(parent->m_rowOffset+ptr(i[idof])->getGlobalId()+ii,
                                                                                          parent->m_colOffset+ptr(j[jdof])->getGlobalId()+jj,
                                                                                          toAssembler(ipos+ii, jpos+jj)));
                            }
                            
                            jpos += jj;
                        }
                    }
                    
                    ipos += ii;
                }
            }
        };*/
        
        AssemblerImplEigenSparseMatrix() : AssemblerBase() {
            m_assembled.resize(1,1);
        }
        ~AssemblerImplEigenSparseMatrix() { }
        
        void init(unsigned int m, unsigned int n) {
            
            //need to set the size of the matrix here which means I need world sizes
            m_assembled.resize(m,n);
            m_tripletList.clear();
        }
        
        void setOffset(unsigned int rowOffset, unsigned int colOffset) {
            AssemblerBase::setOffset(rowOffset, colOffset);
        }
        
        void finalize() {
            
            m_assembled.setFromTriplets(m_tripletList.begin(), m_tripletList.end());
        }
        
        //rewrite --> assemble uses assemble object inline constructor to do all the work.
        template<typename I, typename J, typename Input>
        inline void assemble(I &i, J &j, Input &toAssembler) {
            //std::cout<<"Size: "<<m_assembled.rows()<<" "<<m_assembled.cols()<<"\n";
            //std::cout<<"Not an accepted matrix type to assemble "<<toAssembler<<"\n";
            assembleStruct<I,J,Input>(this, i, j, toAssembler);
            
            //needs to be handled by an object
        }
        
        const  auto & getMatrix() const { return m_assembled; }
        auto & getMatrix() { return m_assembled; }
        
    protected:
        
        //List of triplets
        std::vector<Eigen::Triplet<Precision>> m_tripletList;
        
        //Sparse Matrix
        Eigen::SparseMatrix<Precision, Eigen::RowMajor> m_assembled;
        
    private:
        
    };
   
    class AssemblerImplEigenVector : public AssemblerBase {
        typedef double Precision;
        
    public:
        
        using MatrixType =Eigen::VectorXd;
        using AssemblerBase::m_rowOffset;
        using AssemblerBase::m_colOffset;
        
        template<typename I, typename Input>
        struct assembleStruct {
            
            inline assembleStruct(AssemblerImplEigenSparseMatrix *parent, I i, Input &toAssembler) {
                std::cout<<"Input type not accepted by Assembler \n";
                assert(1==0);
            }
            
        };
        
        template<typename I,  int ROWS>
        struct assembleStruct<I, Eigen::Matrix<Precision, ROWS, 1> > {
            inline assembleStruct(AssemblerImplEigenVector *parent, I &i,
                                  const Eigen::Matrix<Precision, ROWS, 1> &toAssembler) {
                unsigned int localId = 0;
                for(unsigned int idof=0;idof < i.size(); ++idof) {
                        for(unsigned int ii=0; ii<ptr(i[idof])->getNumScalarDOF(); ++ii, ++localId) {
                            parent->m_assembled[parent->m_rowOffset+ptr(i[idof])->getGlobalId()+ii] += toAssembler[localId];
                        }
                    }
                }
        };
        
        AssemblerImplEigenVector() : AssemblerBase() { }
        ~AssemblerImplEigenVector() { }
        
        void init(unsigned int m, unsigned int n) {
            
            //need to set the size of the matrix here which means I need world sizes
            m_assembled.resize(m,1);
            m_assembled.setZero();
        }
        
        void setOffset(unsigned int rowOffset, unsigned int colOffset = 0) {
            AssemblerBase::setOffset(rowOffset, colOffset);
        }
        
        void finalize() {
            
        }
        
        template<typename I, typename Input>
        inline void assemble(I &i, Input &toAssembler) {
            //std::cout<<"Size: "<<m_assembled.rows()<<" "<<m_assembled.cols()<<"\n";
            //std::cout<<"Not an accepted matrix type to assemble "<<toAssembler<<"\n";
            assembleStruct<I, Input>(this, i, toAssembler);
        }
        
        const  Eigen::VectorXd & getVector() const { return m_assembled; }
        Eigen::VectorXd & getMatrix() { return m_assembled; }
        
    protected:
        
        //Sparse Matrix
        Eigen::VectorXd m_assembled;
        
    private:
        
    };
    
    /*deprecated template<>
    inline void AssemblerImplEigenSparseMatrix::assemble<const double>(unsigned int i, unsigned int j, unsigned int size, const double &toAssembler) {
        
        //add double to system as a diagonal matrix
        for(unsigned int ii=0; ii<size; ++ii) {
            m_tripletList.push_back(Eigen::Triplet<Precision>(i+ii,j+ii, toAssembler));
        }
    }*/

    template<typename DataType>
    using AssemblerEigenSparseMatrix = Assembler<DataType, AssemblerImplEigenSparseMatrix>;
    
    template<typename DataType>
    using AssemblerEigenVector = Assembler<DataType, AssemblerImplEigenVector>;
    
    //some utility functions and definitions to make assembling things easier
    //check if something is derived from an assembler
    template<typename T>
    class isAssembler {
    public:
        inline static constexpr bool val() { return false; }
    };

    template<typename DataType, typename Impl>
    class isAssembler<Assembler<DataType, Impl> > {
    public:
        inline static constexpr bool val() { return true; }

    };

    //for matrices
    //equals operator (just selects out for assembler at compile time
    template<typename A, typename B, typename I, typename J, unsigned int Operation=0>
    inline void assign(A &a, B &b, I &&i, J &&j) {
        static_if<isAssembler<typename std::remove_reference<A>::type>::val()>([&](auto f){
            static_if<Operation==0>([&](auto f) { //don't transpose
                f(a).set(i, j, b);
            }).else_([&](auto f) { //transpose
                f(a).set(j,i,b);
            });
        }).else_([&](auto f) {
            f(a) = b;
        });
    }
    
    //for vectors
    //equals operator (just selects out for assembler at compile time
    template<typename A, typename B, typename I>
    inline void assign(A &a, B &b, I &&i) {
        static_if<isAssembler<typename std::remove_reference<A>::type>::val()>([&](auto f){
            f(a).set(i, b);
        }).else_([&](auto f) {
            f(a) = b;
        });
    }

    //for cases where you may need to assign a whole assembled structure to a vector
    template<typename A, typename B>
    inline void assign(A &a, B &b) {
        static_if<isAssembler<typename std::remove_reference<B>::type>::val()>([&](auto f){
            f(a) = (*b);
        }).else_([&](auto f) {
            f(a) = b;
        });
    }
    
}

#endif
