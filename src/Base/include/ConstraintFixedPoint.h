//
//  ConstraintFixedPoint.h
//  Gauss
//
//  Created by David Levin on 4/20/17.
//
//

#ifndef ConstraintFixedPoint_h
#define ConstraintFixedPoint_h

#include <Constraint.h>
#include <DOFParticle.h>
#include <UtilitiesEigen.h>

namespace Gauss {
    template<typename DataType>
    class ConstraintFixedPointImpl
    {
    public:
        
        ConstraintFixedPointImpl(ParticleSystem::DOFParticle<DataType> *q0, Eigen::Vector3d x) {
            m_dofFixed = q0;
            m_p0 = x;
        }
        
        void setFixedPoint(Eigen::VectorXd pos) {
            m_p0 = pos;
        }
        
        ~ConstraintFixedPointImpl() { }
        
        constexpr unsigned int getNumRows() { return 3; }
        //value of constraint (supports vector valued constraint functions for points and what not)
        template<typename Vector>
        inline void getFunction(Vector &f,  const State<DataType> &state, const ConstraintIndex &index) {
            
            Eigen::Vector3d func = m_p0 - mapDOFEigen(*m_dofFixed, state);
            assign(f, func, std::array<ConstraintIndex,1>{{index}});
        }
        
        //get DOFs that this constraint is acting on
        auto & getDOF(unsigned int index) {
            return *m_dofFixed;
        }
        
        //how many DOFs are involved in this constraint
        constexpr unsigned int getNumDOF() const { return 1; }
        
        template <typename Matrix, unsigned int Operation>
        inline void getGradient(Matrix &g,  const State<DataType> &state, const ConstraintIndex &index) {
            Eigen::Matrix3d I;
            I.setIdentity();
            
            assign<Matrix, Eigen::Matrix3d, std::array<ConstraintIndex,1>, std::array<DOFBase<DataType,0>, 1>, Operation>(g, I, std::array<ConstraintIndex,1>{{index}}, std::array<DOFBase<DataType,0>, 1>{{*m_dofFixed}});
        }
        
        
    protected:
        
        Eigen::VectorXd m_p0; //position to fix point at
        ParticleSystem::DOFParticle<DataType> *m_dofFixed; //pointer to the thing I'm fixing in space
        
    private:
    };
    
    template<typename DataType>
    using ConstraintFixedPoint = Constraint<DataType, ConstraintFixedPointImpl<DataType> >;
    
}

#endif /* ConstraintFixedPoint_h */
