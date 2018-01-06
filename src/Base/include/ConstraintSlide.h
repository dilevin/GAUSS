//
//  ConstraintSlide.h
//  Gauss
//
//  Created by David Levin on 1/5/18.
//
//

#ifndef ConstraintSlide_h
#define ConstraintSlide_h


#include <Constraint.h>
#include <DOFParticle.h>
#include <UtilitiesEigen.h>

namespace Gauss {
    template<typename DataType>
    class ConstraintSlideImpl
    {
    public:
        
        ConstraintSlideImpl(ParticleSystem::DOFParticle<DataType> *q0, DataType val, unsigned int dir) {
            
            m_dofFixed = q0;
            m_val = val;
            m_dir = dir;
        }
        
        void setFixedVal(DataType val) {
            m_val = val;
        }
        
        ~ConstraintSlideImpl() { }
        
        constexpr unsigned int getNumRows() { return 1; }
        
        //value of constraint (supports vector valued constraint functions for points and what not)
        template<typename Vector>
        inline void getFunction(Vector &f,  const State<DataType> &state, const ConstraintIndex &index) {
            
            Eigen::Vector3d func = m_val - mapDOFEigen(*m_dofFixed, state)[m_dir];
            assign(f, func, std::array<ConstraintIndex,1>{{index}});
        }
        
        //get DOFs that this constraint is acting on
        auto & getDOF(unsigned int index) {
            return *m_dofFixed;
        }
        
        //how many DOFs are involved in this constraint
        constexpr unsigned int getNumDOF() const { return 1; }
        
        template <typename World, typename Matrix, unsigned int Operation>
        inline void getGradient(Matrix &g,  const World &world, const State<DataType> &state, const ConstraintIndex &index) {
            
            //quick hacky
            static_if<Operation==0>([&](auto f) {
                Eigen::Matrix<DataType,1,3> I;
                I.setZero();
                I(0,m_dir) = 1.0;
                assign<Matrix, Eigen::Matrix<DataType,1,3>, std::array<ConstraintIndex,1>, std::array<DOFBase<DataType,0>, 1>, 0>(g, I, std::array<ConstraintIndex,1>{{index}}, std::array<DOFBase<DataType,0>, 1>{{*m_dofFixed}});
            }).else_([&](auto f) {
                Eigen::Matrix<DataType,3,1> I;
                I.setZero();
                I(m_dir,0) = 1.0;
                assign<Matrix, Eigen::Matrix<DataType,3,1>, std::array<ConstraintIndex,1>, std::array<DOFBase<DataType,0>, 1>, 1>(g, I, std::array<ConstraintIndex,1>{{index}}, std::array<DOFBase<DataType,0>, 1>{{*m_dofFixed}});
            });

            
        }
        
        
    protected:
        
        DataType m_val; //position to fix point at
        unsigned int m_dir;
        ParticleSystem::DOFParticle<DataType> *m_dofFixed; //pointer to the thing I'm fixing in space
        
    private:
    };
    
    template<typename DataType>
    using ConstraintSlide = Constraint<DataType, ConstraintSlideImpl<DataType> >;
    
    //Utility functions to fix a bunch of points
    template<typename World, typename FEMSystem>
    void fixDisplacementDirectionMin(World &world, FEMSystem *system, unsigned int dir = 0, unsigned int dim = 0) {
        //find all vertices with minimum x coordinate and fix DOF associated with them
        auto minX = system->getImpl().getV()(0,dim);
        std::vector<unsigned int> minV;
        
        for(unsigned int ii=0; ii<system->getImpl().getV().rows(); ++ii) {
            
            if(system->getImpl().getV()(ii,dim) < minX) {
                minX = system->getImpl().getV()(ii,dim);
                minV.clear();
                minV.push_back(ii);
            } else if(fabs(system->getImpl().getV()(ii,dim) - minX) < 1e-5) {
                minV.push_back(ii);
            }
        }
        
        //add a bunch of constraints
        for(auto iV : minV) {
            world.addConstraint(new ConstraintSlide<decltype(minX)>(&system->getQ()[iV], 0, dir));
        }
    }

    template<typename World, typename FEMSystem>
    void fixDisplacementDirectionMax(World &world, FEMSystem *system, unsigned int dir = 0, unsigned int dim = 0) {
        //find all vertices with minimum x coordinate and fix DOF associated with them
        auto maxX = system->getImpl().getV()(0,dim);
        std::vector<unsigned int> maxV;
        
        for(unsigned int ii=0; ii<system->getImpl().getV().rows(); ++ii) {
            
            if(system->getImpl().getV()(ii,dim) > maxX) {
                maxX = system->getImpl().getV()(ii,dim);
                maxV.clear();
                maxV.push_back(ii);
            } else if(fabs(system->getImpl().getV()(ii,dim) - maxX) < 1e-5) {
                maxV.push_back(ii);
            }
        }
        
        //add a bunch of constraints
        for(auto iV : maxV) {
            world.addConstraint(new ConstraintSlide<decltype(maxX)>(&system->getQ()[iV], 0, dir));
        }
    }

    
}

#endif /* ConstraintSlide_h */
