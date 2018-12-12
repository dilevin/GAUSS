#include <GaussIncludes.h>
#include <FEMIncludes.h>
#include <Eigen/SparseCholesky>

//Any extra things I need such as constraints
#include <LoubignacIterations.h>
#include <ConstraintFixedPoint.h>

//IGL Viewer
#include <igl/opengl/glfw/viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include <igl/readMESH.h>

//Global variables for UI
igl::opengl::glfw::Viewer viewer;

using namespace Gauss;
using namespace FEM;
using namespace ParticleSystem; //For Force Spring

/* Tetrahedral finite elements */
//typedef scene
typedef PhysicalSystemFEM<double, LinearTet> FEMLinearTets;
typedef World<double, std::tuple<FEMLinearTets *>, std::tuple<ForceSpringFEMParticle<double> *>, std::tuple<ConstraintFixedPoint<double> *> > MyWorld;

//This code from libigl boundary_facets.h
void tetsToTriangles(Eigen::MatrixXi &Fout, Eigen::MatrixXi &T) {
    
    unsigned int simplex_size = 4;
    std::vector<std::vector<int> > allF(
                                        T.rows()*simplex_size,
                                        std::vector<int>(simplex_size-1));
    
    // Gather faces, loop over tets
    for(int i = 0; i< (int)T.rows();i++)
    {
        // get face in correct order
        allF[i*simplex_size+0][0] = T(i,2);
        allF[i*simplex_size+0][1] = T(i,3);
        allF[i*simplex_size+0][2] = T(i,1);
        // get face in correct order
        allF[i*simplex_size+1][0] = T(i,3);
        allF[i*simplex_size+1][1] = T(i,2);
        allF[i*simplex_size+1][2] = T(i,0);
        // get face in correct order
        allF[i*simplex_size+2][0] = T(i,1);
        allF[i*simplex_size+2][1] = T(i,3);
        allF[i*simplex_size+2][2] = T(i,0);
        // get face in correct order
        allF[i*simplex_size+3][0] = T(i,2);
        allF[i*simplex_size+3][1] = T(i,1);
        allF[i*simplex_size+3][2] = T(i,0);
    }
    
    Fout.resize(allF.size(), simplex_size-1);
    for(unsigned int ii=0; ii<allF.size(); ++ii) {
        Fout(ii,0) = allF[ii][0];
        Fout(ii,1) = allF[ii][1];
        Fout(ii,2) = allF[ii][2];
    }
}

int main(int argc, char **argv) {
    std::cout<<"Test Linear FEM \n";
    
    //Setup Physics
    MyWorld world;
    
    //new code -- load tetgen files
    
    Eigen::MatrixXd V;
    Eigen::MatrixXi T;
    
    readTetgen(V, T, dataDir()+"/meshesTetgen/Beam/Beam.node", dataDir()+"/meshesTetgen/Beam/Beam.ele");
    
    Eigen::MatrixXd C;
    Eigen::MatrixXi F;
    Eigen::VectorXd s;
    
    //Eigen::MatrixXd V(4,3);
    //Eigen::MatrixXi T(1,4);
     /*V << 0,0,0,
     1,0,0,
     0,1,0,
     0,0,1;
     
    T << 0, 1, 2, 3;*/
    FEMLinearTets *test = new FEMLinearTets(V,T);
    world.addSystem(test);
    world.finalize();
    
    auto q = mapStateEigen(world);
    q.setZero();
    
    Eigen::VectorXi indices = minVertices(test, 0);
    Eigen::SparseMatrix<double> P = fixedPointProjectionMatrix(indices, *test,world);
    
    AssemblerEigenSparseMatrix<double> K;
    AssemblerEigenVector<double> f;
    
    getStiffnessMatrix(K, world);
    getForceVector(f, world);
    
    Eigen::SparseMatrix<double> Kp = P*(*K)*P.transpose();
    Eigen::VectorXd fp = -P*(*f);
    
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
    solver.compute(Kp);

    if(solver.info()!=Eigen::Success) {
        // decomposition failed
        assert(1 == 0);
        std::cout<<"Decomposition Failed \n";
        exit(1);
    }
    
    //std::cout<<"ANSWER: "<<(P.transpose()*solver.solve(fp)).transpose()<<"\n";
    mapStateEigen<0>(world) = P.transpose()*solver.solve(fp);

    Eigen::MatrixXd stress;
    loubignacIterations(stress, (*K), P, world.getState(), (*f), *test, 1e-7);
    
    s = stress.rowwise().squaredNorm();
    
    tetsToTriangles(F, T);
    
    std::cout<<"Loubignac Iterations: \n"<<s.minCoeff()<<" "<<s.maxCoeff()<<"\n";
    igl::jet(s, s.minCoeff(), 100000000000, C);
    viewer.data().set_mesh(V, F);
    viewer.data().set_colors(C);
    viewer.launch();
    
}

