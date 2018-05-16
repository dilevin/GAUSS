//
//  DOFPair.h
//  Gauss
//
//  Created by David Levin on 5/13/18.
// A pair of DOFS concatanted into a single DOF

#ifndef DOFPair_h
#define DOFPair_h

namespace Gauss {
    
    template<typename DataType, template<typename A, unsigned int B> class DOF1, template<typename C, unsigned int D> class DOF2, unsigned int PropertyIndex=0>
    class DOFPair : public DOFBase<DataType, PropertyIndex> {
        
    public:
        DOFPair() : DOFBase<DataType, PropertyIndex>(0)
        {
            //setup dof localIds
            m_dof1.setIds(0,0);
            m_dof2.setIds(m_dof1.getNumScalarDOF(),m_dof1.getNumScalarDOF());
            
            DOFBase<DataType, PropertyIndex>::m_numScalarDOF = m_dof1.getNumScalarDOF() + m_dof2.getNumScalarDOF();
            
        }
        
        ~DOFPair() { }
        
        inline DOF1<DataType, PropertyIndex> & first() {
            return m_dof1;
        }
        
        inline DOF2<DataType,PropertyIndex> & second() {
            return m_dof2;
        }
        
        inline const DOF1<DataType,PropertyIndex> & first() const {
            return m_dof1;
        }
        
        inline const DOF2<DataType,PropertyIndex> & second() const {
            return m_dof2;
        }
        
        constexpr unsigned int getIndex() { return PropertyIndex; }
        
        inline unsigned int getNumScalarDOF() const {
            return m_dof1.getNumScalarDOF() + m_dof2.getNumScalarDOF();
        }
        
        inline void offsetGlobalId(unsigned int offset) {
            DOFBase<DataType, PropertyIndex>::m_globalId += offset;
            m_dof1.offsetGlobalId(DOFBase<DataType, PropertyIndex>::m_globalId);
            m_dof2.offsetGlobalId(DOFBase<DataType, PropertyIndex>::m_globalId);
        }
        
        inline unsigned int  getNumDOFs() { return 2; }
        
    protected:
        DOF1<DataType,PropertyIndex> m_dof1;
        DOF2<DataType,PropertyIndex> m_dof2;
        
    private:
        
    };
}
#endif /* DOFPair_h */
