#ifndef NEWTON_HPP
#define NEWTON_HPP
#include <vector>

namespace Gauss {
    namespace Optimization {
        
        //this function takes a single Newton Step, it requires 4 functions
        //Energy: returns scalar energy as as function of a Gauss State
        //Gradient: returns the gradient of the energy as  function of a Gauss state
        //Hessian: returns the hessian of the energy as a function of a Gauss state
        //Solve: takes a matrix,A, and a vector,b, as a parameter and returns the solution to Ax = b
        //tol1 is the first order optimality tolerance
        //Note: State is updated inline, so answer is returned in x0.
        //Everything should work with Eigen for now
        template <typename Energy, typename Gradient, typename Hessian, typename Solve, typename PostStepCallback, typename Vector>
        inline bool newtonStep(Energy &f, Gradient &g, Hessian &H, Solve &solver, Vector &x0, PostStepCallback &pscallback, double tol1 = 1e-5) {
            
            //first version is a lazy Newton step with no line search.
            double E0 = f(x0);
            auto p = -solver(H(x0), g(x0));
            double gStep = g(x0).transpose()*p;
            
            //back tracking line search
            double alpha =1;
            double c = 1;
            double rho = 0.5;
            
            while(true) {
                
                pscallback(x0+alpha*p);
                
                if(f(x0+alpha*p) <  E0) {
                    break;
                }
                
                
                alpha = alpha*rho;
            }
            
            x0 = x0+alpha*p;
            
            return (fabs((f(x0) - E0)/E0) < tol1 ? true : false);
            
        }
        
        //poststep callback is a function that takes in x0
        template <typename Energy, typename Gradient, typename Hessian, typename Solve, typename PostStepCallback, typename Vector>
        bool newton(Energy &f, Gradient &g, Hessian &H, Solve &solver, Vector &x0, PostStepCallback &pscallback, double tol1 = 1e-5, unsigned int maxIter = 100) {
            
            
            unsigned int iter = 0;
            
            while(newtonStep(f,g,H, solver, x0, pscallback, tol1) == false && iter < maxIter) {
                iter++;
            }
            
            return true;
        }
        
        //poststep callback is a function that takes in x0
        template<typename Energy, typename Gradient, typename Hessian, typename Solve, typename Vector>
        bool newton(Energy &f, Gradient &g, Hessian &H, Solve &solver, Vector &x0, double tol1 = 1e-5, unsigned int maxIter = 100) {
            
            unsigned int iter = 0;
            
            while(newtonStep(f,g,H, solver, x0, [](auto &x){}, tol1) == false && iter < maxIter) {
                iter++;
            }
            
            return true;
        }
        
    }
    
}
#endif //NEWTON_HPP
