 //
//  ForcePoint.h
//  Gauss
//
//  Created by David Levin on 1/5/18.
//
//

#ifndef ForcePoint_h
#define ForcePoint_h

#include <ForceExternal.h>
#include <DOFParticle.h>
#include <UtilitiesEigen.h>
#include <unsupported/Eigen/AutoDiff>

namespace Gauss {
    namespace ParticleSystem {
        
        template<typename DataType>
        class ForcePointImpl
        {
        public:
            
            ForcePointImpl(DOFBase<DataType> *dof, Eigen::VectorXx<DataType> fPoint) : m_dof(0) {
                m_dof = dof;
                m_fPoint = fPoint;
            }
            
            //forces can give you the energy stored, the force itself and the hessian
            inline DataType getEnergy(State<DataType> &state) {
                
                return -mapDOFEigen(*m_dof, state).dot(m_fPoint);
                
                
            }
            
            //forces always act on at least one DOF of the system this function returns which DOF the are acting on.
            inline auto & getDOF(unsigned int index) {
                
                return *m_dof;
            }
            
            inline unsigned int getNumDOF() { return 1; }
            
            template <typename Vector>
            inline void getForce(Vector &f, State<DataType> &state) {
                
                assign(f, m_fPoint, std::array<DOFBase<DataType,0> , 1>{{*m_dof}});
                
            }
            
            template<typename Vector>
            inline void setForce(Vector &f) {
                m_fPoint = f;
            }
            
            template <typename Matrix>
            inline void getStiffnessMatrix(Matrix &H, State<DataType> &state) {
                
                
            }
            
        protected:
            
            Eigen::VectorXx<DataType> m_fPoint;
            DOFBase<DataType> *m_dof; //pointer to the thing I'm fixing in space
            
        private:
        };
        
        template<typename DataType>
        using ForcePoint = Force<DataType, ForcePointImpl<DataType> >;
        
    }
}

#endif /* ForcePoint_h */
