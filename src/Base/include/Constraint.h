#ifndef _CONSTRAINT_H
#define _CONSTRAINT_H

namespace Gauss {
    
    //check for update method
    template <typename T, typename World>
    inline auto callUpdate (T &t, World &world, int i) -> decltype( t.update(world) )
    { t.update(world); }
    
    template <typename T, typename World>
    inline auto callUpdate (T &t, World &world, long i)
    {  }
    
    class ConstraintIndex {
    public:
        
        inline ConstraintIndex(unsigned int localId, unsigned int globalId, unsigned int numRows) {
            setIds(localId, globalId, numRows);
        }
        
        inline void setIds(unsigned int localId, unsigned int globalId, unsigned int numRows) {
            m_localId = localId;
            m_globalId = globalId;
            m_numRows = numRows;
        }
        
        inline void setNumRows(unsigned int numRows) { m_numRows = numRows; }
        
        inline unsigned int getNumScalarDOF() {
            return m_numRows;
        }
        
        inline unsigned int getNumScalarDOF() const {
            return m_numRows;
        }
        
        inline unsigned int getLocalId() const {
            return m_localId;
        }
        
        inline unsigned int getGlobalId() const {
            return m_globalId;
        }
        
        inline void setGlobalId(unsigned int globalId) {
            m_globalId = globalId;
        }
        
        virtual inline void offsetGlobalId(unsigned int offset) {
            m_globalId += offset;
        }

        
        
    protected:
        
        unsigned int m_localId;
        unsigned int m_globalId;
        unsigned int m_numRows;
        
    private:
    };
    
    //Constraints are in the form f(q(t)) = b(t) where b is the constant term in the constraint function.
    //q are the degrees of freedom of the system
    //q and b can be parameterized by time.
    //Constraints need to be able to return f, b, dF/dq, db/dt and the constraint error (f(q(t)) - b(t)).
    template<typename DataType, typename Impl >
    class Constraint
    {
    public:
        template<typename ...Params>
        Constraint(Params ... params) : m_impl(params...), m_index(0,0,m_impl.getNumRows()) {}
        ~Constraint() { }
        
        //value of constraint (supports vector valued constraint functions for points and what not)
        template<typename Vector>
        inline void getFunction(Vector &f,  const State<DataType> &state) {
            m_impl.getFunction(f, state, m_index);
        }
        
        inline ConstraintIndex & getIndex() { return m_index; }
        
        inline unsigned int getNumRows() { m_index.setNumRows(m_impl.getNumRows()); return m_impl.getNumRows(); }
            
        //get DOFs that this constraint is acting on
        auto & getDOF(unsigned int index) {
            return m_impl.getDOF(index);
        }
        
        //how many DOFs are involved in this constraint
        unsigned int getNumDOF() { return m_impl.getNumDOF(); }
        
        template <typename World, typename Matrix>
        inline void getError(Matrix &g,  World &world, const State<DataType> &state) {
            m_impl.getError(g, state, m_index);
        }
        
        template <typename World, typename Vector>
        inline void getB(Vector &b,  World &world, const State<DataType> &state) {
            m_impl.getB(b, state, m_index);
        }
        
        template <typename World, typename Vector>
        inline void getDbDt(Vector &g,  World &world, const State<DataType> &state) {
            m_impl.template getDbDt<World, Matrix,Operation>(g, world, state, m_index);
        }
        
        template <typename World, typename Matrix, unsigned int Operation=0>
        inline void getGradient(Matrix &g,  World &world, const State<DataType> &state) {
            m_impl.template getGradient<World, Matrix,Operation>(g, world, state, m_index);
        }
        
        
        template<typename World>
        inline void update(World &world) {
            callUpdate(m_impl, world, 0);
        }
        
        inline auto & getImpl() { return m_impl; }
    protected:
        
        Impl m_impl;
        ConstraintIndex m_index;
    
    private:
    };

}

#endif
