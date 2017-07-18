//
//  AssemblerMVP.h
//  Gauss
//
//  Created by David Levin on 7/4/17.
//
//

#ifndef AssemblerMVP_h
#define AssemblerMVP_h

//An assembler that performs a global matrix vector product without assembling the full global matrix
#include <Assembler.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <UtilitiesEigen.h>

namespace Gauss {
    class AssemblerMVPImplEigen : public AssemblerBase {
        typedef double Precision;
        
    public:
        
        using MatrixType = Eigen::VectorXx<Precision>;
        
        using AssemblerBase::m_rowOffset;
        using AssemblerBase::m_colOffset;
        
        template<typename I, typename J, typename Input>
        struct assembleStruct {
            
            inline assembleStruct(AssemblerMVPImplEigen *parent, I i, J j, Input &toAssembler) {
                std::cout<<"Input type not accepted by Assembler \n";
                assert(1==0);
            }
            
            inline assembleStruct(AssemblerMVPImplEigen *parent, I i, Input &toAssembler) {
                std::cout<<"Input type not accepted by Assembler \n";
                assert(1==0);
            }
        };
        
        template<typename I, typename J>
        struct assembleStruct<I, J, const double> {
            inline assembleStruct(AssemblerMVPImplEigen *parent, I &i, J &j, const double &toAssembler) {
                
                //add double to system as a diagonal matrix
                for(unsigned int idof=0;idof < i.size(); ++idof) {
                    for(unsigned int jdof=0;jdof < j.size(); ++jdof) {
                        for(unsigned int ii=0; ii<ptr(i[idof])->getNumScalarDOF(); ++ii) {
                            
                            //do the matrix vector product and add it to mvp
                            parent->m_b[parent->m_rowOffset+ptr(i[idof])->getGlobalId()+ii] += toAssembler*(*parent->m_x)[parent->m_colOffset+ptr(j[jdof])->getGlobalId()+ii];
                            //parent->m_tripletList.push_back(Eigen::Triplet<Precision>(parent->m_rowOffset+ptr(i[idof])->getGlobalId()+ii,
                                                                                      //parent->m_colOffset+ptr(j[jdof])->getGlobalId()+ii,
                                                                                      //toAssembler));
                        }
                    }
                }
            }
        };
        
        template<typename I, typename J, int ROWS, int COLS>
        struct assembleStruct<I, J, Eigen::Matrix<double, ROWS, COLS> > {
            inline assembleStruct(AssemblerMVPImplEigen *parent, I &i, J &j, const Eigen::Matrix<double, ROWS, COLS> &toAssembler) {
                
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
                                //parent->m_tripletList.push_back(Eigen::Triplet<Precision>(parent->m_rowOffset+ptr(i[idof])->getGlobalId()+ii,
                                  //                                                        parent->m_colOffset+ptr(j[jdof])->getGlobalId()+jj,
                                  //                                                        toAssembler(ipos+ii, jpos+jj)));
                                parent->m_b[parent->m_rowOffset+ptr(i[idof])->getGlobalId()+ii] += toAssembler(ipos+ii, jpos+jj)*(*parent->m_x)[parent->m_colOffset+ptr(j[jdof])->getGlobalId()+jj];
                            }
                            
                            jpos += jj;
                        }
                    }
                    
                    ipos += ii;
                }
            }
        };
        
        AssemblerMVPImplEigen() : AssemblerBase() {
            
        }
        
        ~AssemblerMVPImplEigen() { }
        
        void init(unsigned int m, unsigned int n) {
            
            m_n = n;
            m_b.resize(m,1);
            m_b.setZero();
        }
        
        //for matrix vector product assembler
        void setX(MatrixType &x) {
            
            if(x.rows() != m_n) {
                std::cout<<"Dimensions of x do not match number of columns in matrix \n";
                exit(0);
            }
            
            m_x = &x;
        }
        
        void setOffset(unsigned int rowOffset, unsigned int colOffset) {
            AssemblerBase::setOffset(rowOffset, colOffset);
        }
        
        void finalize() {
            
        }
        
        //rewrite --> assemble uses assemble object inline constructor to do all the work.
        template<typename I, typename J, typename Input>
        inline void assemble(I &i, J &j, Input &toAssembler) {
            //std::cout<<"Size: "<<m_assembled.rows()<<" "<<m_assembled.cols()<<"\n";
            //std::cout<<"Not an accepted matrix type to assemble "<<toAssembler<<"\n";
            assembleStruct<I,J,Input>(this, i, j, toAssembler);
            
            //needs to be handled by an object
        }
        
        const  auto & getMatrix() const { return m_b; }
        auto & getMatrix() { return m_b; }
        
        template<typename VectorType>
        inline AssemblerMVPImplEigen & operator*=(VectorType &x) { setX(x); return *this; }
        
    protected:
    
        MatrixType m_b; //the matrix vector product
        MatrixType *m_x; //pointer to vector to multiply
        unsigned int m_n;
        
    private:
        
    };
    
    template<typename DataType>
    using AssemblerMVPEigen = Assembler<DataType, AssemblerMVPImplEigen>;


}
#endif /* AssemblerMVP_h */
