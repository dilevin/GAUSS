#include <UtilitiesEigen.h>

#include <FEMIncludes.h>
#include <DOFStress.h>
#include <EnergyLoubignac.h>
#include <ConstraintFixedPoint.h>

using namespace Gauss;
using namespace FEM;
using namespace ParticleSystem; //For Force Spring


namespace Gauss {
    namespace FEM {

        template<typename DataType>
        using LoubignacTet = Element<DataType, 4, QuadratureTetLinear, EnergyKineticNonLumped, EnergyLoubignac, BodyForceGravity, ShapeFunctionLinearTet>;
        
        typedef PhysicalSystemFEM<double, LoubignacTet> FEMLoubignacTet;
        
        template<typename DataType>
        inline void getCauchyStressVector(Eigen::MatrixXx<DataType> &stresses, const State<DataType> &u, const PhysicalSystemFEM<DataType, LinearTet> &fem) {
            
            stresses.resize(fem.getGeometry().second.rows(), 6);
            
            Eigen::Matrix33x<DataType> stress;
            
            for(unsigned int ii=0; ii<fem.getGeometry().second.rows(); ++ii) {
                fem.getImpl().getElement(ii)->getCauchyStress(stress, Vec3d(0,0,0), u);
                stresses(ii,0) = stress(0,0);
                stresses(ii,1) = stress(1,1);
                stresses(ii,2) = stress(2,2);
                stresses(ii,3) = stress(1,2);
                stresses(ii,4) = stress(0,2);
                stresses(ii,5) = stress(0,1);
            }
        
            
        }
        
        template<typename DataType,
                 typename DerivedElementStress,
                 typename DerivedVertStress>
        inline void averageStressesOntoVertices(Eigen::MatrixBase<DerivedVertStress> &sv,
                                                const Eigen::MatrixBase<DerivedElementStress> &se,  const PhysicalSystemFEM<DataType, LinearTet> &fem) {
            
            sv.setZero();
            
            auto T = fem.getImpl().getGeometry().second;
            auto V = fem.getImpl().getGeometry().first;
            
            Eigen::VectorXi N;
            N.resize(V.rows());
            N.setZero();
            
            for(unsigned int ii=0; ii< T.rows(); ii++) {
                
                //only for Tets
                for(unsigned int jj=0; jj<4; ++jj) {
                    N[T(ii,jj)] += 1; //get normalization values for each vertex
                }
            }
            
            for(unsigned int ii=0; ii< T.rows(); ii++) {
                
                for(unsigned int jj=0; jj<4; ++jj) {
                    
                    //average and convert format here (you can do this in one line of eigen code)
                    sv(6*T(ii,jj)) +=   se(ii, 0)/static_cast<double>(N(T(ii,jj)));
                    sv(6*T(ii,jj)+1) += se(ii, 1)/static_cast<double>(N(T(ii,jj)));
                    sv(6*T(ii,jj)+2) += se(ii, 2)/static_cast<double>(N(T(ii,jj)));
                    sv(6*T(ii,jj)+3) += se(ii, 3)/static_cast<double>(N(T(ii,jj)));
                    sv(6*T(ii,jj)+4) += se(ii, 4)/static_cast<double>(N(T(ii,jj)));
                    sv(6*T(ii,jj)+5) += se(ii, 5)/static_cast<double>(N(T(ii,jj)));
                }
            }
            
        }
        
        //Loubignac Iterations for generating smooth stress fields from FEM results
        template<typename DataType,
                typename DerivedK,
                typename DerivedP,
                typename DerivedF>
        inline void loubignacIterations(Eigen::MatrixXx<DataType> &vertexStress,
                    const Eigen::SparseMatrixBase<DerivedK> &K, const Eigen::SparseMatrix<DerivedP> &P,
                    const State<DataType> &u,
                    const Eigen::MatrixBase<DerivedF> &f,
                    const PhysicalSystemFEM<DataType, LinearTet> &fem, DataType tol,
                    unsigned int numIterations = 2000) {
            
            
            typedef World<double, std::tuple<FEMLoubignacTet *>, std::tuple<ForceSpringFEMParticle<double> *>, std::tuple<ConstraintFixedPoint<double> *> > LWorld;
            
            //Build FE system to help integrate stresses
            std::cout<<"Loubignac Iterations\n";
            
            
            State<DataType> smoothStresses(fem.getGeometry().first.rows()*6, fem.getGeometry().first.rows()*120);
            
            FEMLoubignacTet *test = new FEMLoubignacTet(fem.getGeometry().first, fem.getGeometry().second);
            LWorld world;
            world.addSystem(test);
            world.finalize();
            
            
            //Use eigen solvers for now because I'm lazy and I haven't quite thought about how I want to handle variable solvers
            Eigen::SparseMatrix<DataType> Kp = P*K*P.transpose();
            Eigen::VectorXx<DataType> fp = P*f;
            
            Eigen::SimplicialLDLT<Eigen::SparseMatrix<DataType> > solver;
            solver.compute(Kp);
            
            Eigen::VectorXx<DataType> ps = 0*fp;
            Eigen::VectorXx<DataType> du = 0*f;
            Eigen::MatrixXx<DataType> elementStress;
            
            auto smoothStressVector = mapStateEigen<0>(smoothStresses);
    
            AssemblerParallel<DataType, AssemblerEigenVector<DataType> > stressAssembler;
            
            //Iterations
            for(unsigned int ii=0; ii<numIterations; ++ii) {
                
                //Compute per element Stresses
                getCauchyStressVector(elementStress, u, fem);
                
                //Average stresses onto vertices
                averageStressesOntoVertices(smoothStressVector, elementStress, fem);
                
                //Integrate smooth stresses
                getInternalForceVector(stressAssembler, *test, smoothStresses);
               
                //Solve for du
                ps = -fp - P*(*stressAssembler);
                
                du = P.transpose()*solver.solve(ps);
                
                mapStateEigen<0>(u) += du;
                
                if(du.norm() < tol) {
                    std::cout<<"Terminated with norm below "<<tol<<" after "<<ii<<" iterations \n";
                    break;
                }
                
               
            }
            
            std::cout<<"Terminated with norm  "<<du.norm()<<"\n";
            vertexStress = (Eigen::Map<Eigen::MatrixXx<DataType> >(smoothStressVector.data(), 6, fem.getGeometry().first.rows())).transpose();
            
        }
    }
}
