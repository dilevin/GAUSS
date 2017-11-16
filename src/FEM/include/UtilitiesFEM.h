//
//  UtilitiesFEM.h
//  Gauss
//
//  Created by David Levin on 8/23/17.
//
//

#ifndef UtilitiesFEM_h
#define UtilitiesFEM_h

#include <GaussIncludes.h>
#include <FEMIncludes.h>
#include <ForceSpring.h>
#include <UtilitiesEigen.h>
#include <complex>

//A collection of useful FEM utilities

namespace Gauss {
    namespace FEM {
        #ifdef GAUSS_SPECTRA
        //linear modal analysis
        //returns a pair wherein the first value is a matrix of eigenvectors and the second is vector of corresponding vibrational frequencies
        template<typename World>
        auto linearModalAnalysis(World &world, unsigned int numModes) {
            
            //build mass and stiffness matrices
            AssemblerParallel<double, AssemblerEigenSparseMatrix<double> > mass;
            AssemblerParallel<double, AssemblerEigenSparseMatrix<double> > stiffness;
            
            getMassMatrix(mass, world);
            getStiffnessMatrix(stiffness, world);
            
            auto eigs = generalizedEigenvalueProblem((*stiffness), (*mass), numModes);
            
            //convert to vibrational frequencies
            //hack because generalised eignvalue solver only returns real values
            for(unsigned int ii=0; ii<eigs.second.rows(); ++ii) {
                eigs.second[ii] = std::sqrt(std::fabs(eigs.second[ii]));
            }
            
            return eigs;
        }
        #endif
        
        //functor for getting position of a DOF
        template <typename DataType, typename DOF>
        class PositionFEMEigen {
        public:
            inline PositionFEMEigen(DOF *dof = nullptr, unsigned int vId = 0, Eigen::MatrixXd *V = 0) {
                m_dof = dof;
                m_V  = V;
                m_vId = vId;
            }
            inline Eigen::Vector3x<DataType> operator()(const State<DataType> &state) const {
                return mapDOFEigen(*m_dof, state) + m_V->row(m_vId).transpose();
            }
            
            inline DOF * getDOF() { return m_dof; }
        protected:
            DOF *m_dof;
            Eigen::MatrixXd *m_V;
            unsigned int m_vId;
        };
        
        //FEM-PARTICLE SPRING
        template<typename DataType>
        using PosFEM = PositionFEMEigen<DataType, ParticleSystem::DOFParticle<DataType> >;
        
        template<typename DataType>
        using ForceSpringFEMParticle = Force<DataType, ParticleSystem::ForceSpringImpl<DataType, PosFEM<DataType>, ParticleSystem::PosParticle<DataType> > >;
    }
}

#endif /* UtilitiesFEM_h */
