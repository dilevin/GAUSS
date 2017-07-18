//
//  SolverCG.h
//  Gauss
//
//  Created by David Levin on 6/28/17.
//
//

#ifndef SolverCG_h
#define SolverCG_h

namespace Gauss {
  
    template<typename DataType, typename VectorType>
    class SolverCG
    {
    public:
        SolverCG(DataType rtol=1e-8) { m_rtol = rtol; }
        
        //Matrix-Vector-Product Operation, RHS
        template<typename MVPFunc, typename PCFunc = const std::function<VectorType&(VectorType &)> >
        void solve(VectorType &x, MVPFunc mvp, VectorType &rhs, unsigned int maxIter = 1e8, PCFunc &pc = [](VectorType &x) -> VectorType & {return x;});
        
    protected:
        
        DataType m_rtol;
        VectorType m_r[2], m_p[2], m_z[2];
        
    private:
    };
    
    //Replaces initial guess with solution
    template<typename DataType, typename VectorType>
    template<typename MVPFunc, typename PCFunc>
    void SolverCG<DataType, VectorType>::solve(VectorType &x, MVPFunc mvp, VectorType &rhs, unsigned int maxIter, PCFunc &pc) {
        
        //CG Loop
        m_r[0].resize(x.rows(),1);
        m_p[0].resize(x.rows(),1);
        m_z[0].resize(x.rows(),1);
        
        
        //use x as the initial guess
        
        m_r[0] = rhs - mvp(x);
        m_z[0] = pc(m_r[0]);
        m_p[0] = m_z[0];
        
        DataType alpha, beta;
        
        unsigned int k=0;
        unsigned int kp1 = 0;
        for(unsigned int ii=0; ii<maxIter; ++ii) {
            kp1 = 1-k;
            alpha = m_r[k].dot(m_z[k])/(mvp(m_p[k]).dot(m_p[k]));
            x += alpha*m_p[k];
            m_r[kp1] = m_r[k] - alpha*mvp(m_p[k]);
            
            if(m_r[kp1].norm() < m_rtol)
                break;
            
            m_z[kp1] = pc(m_r[kp1]);
            beta = m_z[kp1].dot(m_r[kp1])/(m_z[k].dot(m_r[k]));
            
            
            m_p[kp1] = m_z[kp1] + beta*m_p[k];
            
            k = 1-k; //swap buffer
            
        }
        
    }
}
#endif /* SolverCG_h */
