//
//  DOFList.h
//  Gauss
//
//  Created by David Levin on 3/21/17.
//
//

#ifndef DOFList_h
#define DOFList_h

#include  <vector>
namespace Gauss {
    namespace FEM {
        template<typename DataType, template<typename A, unsigned int B> class DOF, unsigned int PropertyIndex = 0>
        class DOFList : public DOFBase<DataType, PropertyIndex> {
        
        public:
            DOFList(unsigned int numDOFs) : DOFBase<DataType, PropertyIndex>(0), m_dofList(numDOFs)
            {
                //setup dof localIds
                unsigned int localId = 0;
                for(auto &dof : m_dofList) {
                    dof.setIds(localId, localId);
                    localId += dof.getNumScalarDOF();
                }
            }
            
            ~DOFList() { }
            
            inline DOF<DataType, PropertyIndex> & operator[](unsigned int index) {
                assert(index < m_dofList.size());
                return m_dofList[index];
            }
            
            constexpr unsigned int getIndex() { return PropertyIndex; }
            
            
            //inline std::pair<DataType *, unsigned int> getPtr(const State<DataType> &state) {
             //   return state.template getStatePtr<PropertyIndex>(DOFBase<DataType, PropertyIndex>::m_globalId);
            //}
            
            inline unsigned int getNumScalarDOF() const {
                return m_dofList[0].getNumScalarDOF()*m_dofList.size();
            }
            
            inline void offsetGlobalId(unsigned int offset) {
                DOFBase<DataType, PropertyIndex>::m_globalId += offset;
                
                //update all the dofs in the dof list
                for(auto &dof : m_dofList) {
                    dof.offsetGlobalId(DOFBase<DataType, PropertyIndex>::m_globalId);
                }
            }
            
            inline unsigned int &getNumDOFs() { return m_dofList.size(); }

        protected:
            std::vector<DOF<DataType, PropertyIndex> > m_dofList;
        private:
        };
    }
}
#endif /* DOFList_h */
