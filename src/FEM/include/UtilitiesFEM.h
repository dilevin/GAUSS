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
#include <Constraint.h>

//A collection of useful FEM utilities

namespace Gauss {
    namespace FEM {
        
        //linear modal analysis
        //returns a pair wherein the first value is a matrix of eigenvectors and the second is vector of corresponding vibrational frequencies
        template<typename World>
        auto linearModalAnalysis(World &world, unsigned int numModes) {
            
            //build mass and stiffness matrices
            AssemblerParallel<double, AssemblerEigenSparseMatrix<double> > mass;
            AssemblerParallel<double, AssemblerEigenSparseMatrix<double> > stiffness;
            
            getMassMatrix(mass, world);
            getStiffnessMatrix(stiffness, world);
            
            auto eigs = generalizedEigenvalueProblem((*stiffness), (*mass), numModes,   1e-6);
            
            //convert to vibrational frequencies
            //hack because generalised eignvalue solver only returns real values
            for(unsigned int ii=0; ii<eigs.second.rows(); ++ii) {
                eigs.second[ii] = std::sqrt(std::fabs(eigs.second[ii]));
            }
            
            return eigs;
        }
        
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
        using PosFEM = PositionFEMEigen<DataType, DOFParticle<DataType> >;
        
        template<typename DataType>
        using ForceSpringFEMParticle = Force<DataType, ParticleSystem::ForceSpringImpl<DataType, PosFEM<DataType>, ParticleSystem::PosParticle<DataType> > >;
        
        //some convenience functions for getting various, useful, assembled matrices
        //FEM is assumed to be an FEM system
        //x is a nx3 matrix of points in space (n is # of points), we're going to build the shape function matrix evaluated at each point
        //element[i] is a n-vector that stores the index of the element containing the ith vertex in the embedded mesh 
        template<typename Matrix, typename Vector, typename FEM>
        void getShapeFunctionMatrix(Matrix &N, Vector &element, Eigen::MatrixXd &x, FEM &fem) {
            
            double y[3];
            ConstraintIndex cIndex(0, 0, 3); //ConstraintIndices help the assembler understand how rows are handled globally
            element.resize(x.rows(),1);
                
            ASSEMBLEMATINIT(N, 3*x.rows(), fem.getQ().getNumScalarDOF());
//            int vert_count = 0;
            bool counted = false;
            for(unsigned int ii=0; ii<x.rows(); ++ii) {
                
                // to find the closest element, if the vert is not inside any
                counted = false;
                double smallest_bary = -1;
                int smallest_bary_element = 0;
                double y_smallest_bary[3];

                for(unsigned int jj=0; jj<fem.getElements().size(); ++jj) {
                    //compute shape function matrix for an element, if it has negative values we're outside so carry on.
                    y[0]= x(ii,0);
                    y[1]= x(ii,1);
                    y[2]= x(ii,2);
                    
                    auto Jmat = fem.getElements()[jj]->N(y);
//                    std::cout<<Jmat.cols()<<std::endl;
//                    std::cout<<Jmat.rows()<<std::endl;
//                    std::cout<<Jmat<<std::endl;
//                    if (Jmat.minCoeff() >= 0 && Jmat.maxCoeff() > 0) {
                    
                    if (Jmat.minCoeff() >= -1e-12) {
                            //use the assembler to add it to the current matrix
                        N.set(std::array<ConstraintIndex,1>{{cIndex}}, fem.getElements()[jj]->q(), Jmat);
                        element[ii] = jj;
//                        ++vert_count;
                        counted = true;
                        break;
                    }
                    else
                    {
                        if(Jmat.minCoeff() > smallest_bary)
                        {
                            smallest_bary = Jmat.minCoeff();
//                            y_smallest_bary = y;
                            y_smallest_bary[0]= y[0];
                            y_smallest_bary[1]= y[1];
                            y_smallest_bary[2]= y[2];
                            
//                            N.set(std::array<ConstraintIndex,1>{{cIndex}}, fem.getElements()[jj]->q(), Jmat);
                            element[ii] = jj;
                        }
                    }
                    
                    
                }
                
                
                if(counted == false)
                {
                    auto Jmat = fem.getElements()[element[ii]]->N(y_smallest_bary);
                    N.set(std::array<ConstraintIndex,1>{{cIndex}}, fem.getElements()[element[ii]]->q(), Jmat);
//                    std::cout<<y<<std::endl;
//                    std::cout<<ii<<std::endl;
                }
                cIndex.offsetGlobalId(3); //increment global id by 3 for the next point
            }
//            std::cout<<x.rows()<<std::endl;
//            std::cout<<vert_count<<std::endl;;
//            std::cout<<fem.getElements().size();
            
            ASSEMBLEEND(N);
        }
    }
}

#endif /* UtilitiesFEM_h */
