#ifndef NEWTON_HPP
#define NEWTON_HPP
#include <vector>
#include <climits>

namespace Gauss {
    namespace Optimization {
        
        //Implementation of backtracking linesearch. Parameter values chosen based on Numerical Optimization by Nocedal and Wright
        //Inputs:
        //  Energy - a function/functor of the form f(x) that returns the energy associated with the  iterate x
        //  Hessian - a function/functor of the form H(x) that returns the hessian of f(x) associated with the iterate x
        //  ConstraintEq - a function/functor of the form C(x) that returns the vector valued constraint function associated with iterate x
        // JacobianEq -   a function/functor of the form J(x) that returns the gradient of C(x) associated with the iterate x
        // Direction - a function/functor of the form solveR(H, g, Aeq, ceq) that produces a direction for the current line search
        // x0 - initial point for the solver. For constrained solves, x0 is the stacked vector of primal and dual variables. Must be of size #primal + #constraints
        // scallback - a callback function/functor of the form s(x) that is called every time the iterate, x, is updated.
        template <typename Energy, typename Gradient, typename Hessian,
                  typename ConstraintEq, typename JacobianEq,
                  typename Direction, typename StepCallback, typename Vector>
        inline bool backTrackingLinesearch(Vector &x0, Energy &f, Gradient &g, Hessian &H, ConstraintEq &ceq, JacobianEq &Aeq,
                                 Direction &solver, StepCallback &scallback, double tol1 = 1e-5) {
            
            scallback(x0);
            
            Vector gradient;
            auto J = *Aeq(x0);
            auto b = *ceq(x0);
            assign(gradient, g(x0));
            
            double E0 = f(x0) + 100000.0*(J*x0.head(gradient.rows()) - b).norm(); //merit function
            Vector p = solver(H(x0), gradient, ceq(x0), Aeq(x0), x0);
            double gStep = gradient.transpose()*p.head(gradient.rows());
            
            //back tracking line search
            double alpha = 1;
            double c = 1e-8; //from Nocedal and Wright pg 31
            double rho = 0.5;
            
            //for(unsigned int ii=0; ii < 20; ++ii) {
            while(true) {
                
                scallback(x0+alpha*p);
                
                if((f(x0+alpha*p)+ 100000.0*(J*(x0+alpha*p).head(gradient.rows()) - b).norm() - E0) < tol1) {
                    break;
                }
                
                std::cout<<f(x0+alpha*p)+ 100000.0*(J*(x0+alpha*p).head(gradient.rows()) - b).norm()<<" "<<E0 + c*alpha*gStep<<" "<<E0<<" "<<alpha<<"\n";
                alpha = alpha*rho;
                
                if(alpha < 1e-8) {
                  alpha = 0;
                }
            }
            
            x0 = x0+alpha*p;
            
            return ((fabs(f(x0) - E0)) < tol1  ? true : false);
            
        }
        
        //Outer loop for a linesearch based optimization
        //  Energy - a function/functor of the form f(x) that returns the energy associated with the  iterate x
        //  ConstraintEq - a function/functor of the form C(x) that returns the vector valued constraint function associated with iterate x
        // JacobianEq -   a function/functor of the form J(x) that returns the gradient of C(x) associated with the iterate x
        // LineSearch - A function/function of the form  linesearch(f,g,H, ceq, Aeq, x0, tol1) that computes a search direction then does linesearch
        // x0 - initial point for the solver. For constrained solves, x0 is the stacked vector of primal and dual variables. Must be of size #primal + #constraints
        //tol1 - convergence tolerane
        template <typename Energy, typename Gradient, typename Hessian,
        typename ConstraintEq, typename JacobianEq, typename Vector, typename Linesearch>
        bool optimizeWithLineSearch(Vector &x0, Energy &f, Gradient &g, Hessian &H, ConstraintEq &ceq, JacobianEq &Aeq, Linesearch &linesearch, double tol1 = 1e-5, unsigned int maxIter = UINT_MAX) {
            
            unsigned int iter = 0;
            
            while(iter < maxIter && linesearch(x0, f,g,H, ceq, Aeq, tol1) == false) {
                iter++;
            }
            
            return true;
            
        }
        
        /*
         backtracking linesearch(Energy, Gradient, Hesssian, Aeq, beq, Aineq, bineq, Solve, , callback)
         
         //
         newton(world, assemblers) do Newton's seaech
         {
            lineSearch(for()
         }
         gradient(assember, world)
         //newtonStep(gauss energy, force, stiffness, constraint
        Define a bunch of things that equate to a newton's solve
        //results Struct = newtonStep(function, function, function, function, vector, function, tolerance
         
         */
        //this function takes a single Newton Step, it requires 4 functions
        //Energy: returns scalar energy as as function of a Gauss State
        //Gradient: returns the gradient of the energy as  function of a Gauss state
        //Hessian: returns the hessian of the energy as a function of a Gauss state
        //Solve: takes a matrix,A, and a vector,b, as a parameter and returns the solution to Ax = b
        //tol1 terminate when relative change in function vale falls below this tolerance
        //Note: State is updated inline, so answer is returned in x0.
        //Everything should work with Eigen for now
        template <typename Energy, typename Gradient, typename Hessian, typename Solve, typename PostStepCallback, typename Vector>
        inline bool newtonStep(Energy &f, Gradient &g, Hessian &H, Solve &solver, Vector &x0, PostStepCallback &pscallback, double tol1 = 1e-5) {
            
             //std::cout<<"NORM2: "<<x0.norm()<<"\n";
            //first version is a lazy Newton step with no line search.
            
            pscallback(x0);
            double E0 = f(x0);
            Vector p = -solver(H(x0), g(x0));
            double gStep = g(x0).transpose()*p;
            
            //back tracking line search
            double alpha =1;
            double c = 1e-4; //from Nocedal and Wright pg 31
            double rho = 0.5;
            
            while(true) {
                
                pscallback(x0+alpha*p);
                
                if((f(x0+alpha*p) - E0 - c*alpha*gStep) < tol1) {
                    break;
                }
                
                //std::cout<<f(x0+alpha*p)<<" "<<E0 + c*alpha*gStep<<" "<<E0<<" "<<alpha<<"\n";
                alpha = alpha*rho;
                
                //if(alpha < 1e-8) {
                  //  alpha = 0;
                //}
            }
            
            x0 = x0+alpha*p;
            
            return (fabs(f(x0) - E0) < tol1 ? true : false);
            
        }
        
        //poststep callback is a function that takes in x0
        template <typename Energy, typename Gradient, typename Hessian, typename Solve, typename PostStepCallback, typename Vector>
        bool newton(Energy &f, Gradient &g, Hessian &H, Solve &solver, Vector &x0, PostStepCallback &pscallback, double tol1 = 1e-5, unsigned int maxIter = 100) {
            
            
            unsigned int iter = 0;
            
            while(iter < maxIter && newtonStep(f,g,H, solver, x0, pscallback, tol1) == false) {
                iter++;
            }
            
            return true;
        }
        
        //poststep callback is a function that takes in x0
        template<typename Energy, typename Gradient, typename Hessian, typename Solve, typename Vector>
        bool newton(Energy &f, Gradient &g, Hessian &H, Solve &solver, Vector &x0, double tol1 = 1e-5, unsigned int maxIter = 100) {
            
            unsigned int iter = 0;
            
            while(iter < maxIter && newtonStep(f,g,H, solver, x0, [](auto &x){}, tol1) == false) {
                iter++;
            }
            
            return true;
        }
        
    }
    
}
#endif //NEWTON_HPP
