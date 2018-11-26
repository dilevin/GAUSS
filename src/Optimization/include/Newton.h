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
        // p is the direction along which to perform line search
        // x0 - initial point for the solver. For constrained solves, x0 is the stacked vector of primal and dual variables. Must be of size #primal + #constraints
        // scallback - a callback function/functor of the form s(x) that is called every time the iterate, x, is updated.
        template <typename Energy, typename Gradient, typename Hessian,
                  typename ConstraintEq, typename JacobianEq,
                  typename StepCallback, typename Vector>
        inline bool backTrackingLinesearch(Vector &x0, Vector &p, Energy &f, Gradient &g, Hessian &H, ConstraintEq &ceq, JacobianEq &Aeq, StepCallback &scallback, double tol1 = 1e-5) {
            
            scallback(x0);
            
            Vector gradient;
            auto J = *Aeq(x0);
            auto b = *ceq(x0);
            assign(gradient, g(x0));
            
            double E0 = f(x0) + 100000.0*(J*x0.head(gradient.rows()) - b).norm(); //merit function
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
                
                //std::cout<<f(x0+alpha*p)+ 100000.0*(J*(x0+alpha*p).head(gradient.rows()) - b).norm()<<" "<<E0 + c*alpha*gStep<<" "<<E0<<" "<<alpha<<"\n";
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
        // Direction solver - function which computes search direction from energy, gradient, hessian  and consrtaints
        // LineSearch - A function/function of the form  linesearch(f,g,H, ceq, Aeq, x0, tol1) that computes a search direction then does linesearch
        // x0 - initial point for the solver. For constrained solves, x0 is the stacked vector of primal and dual variables. Must be of size #primal + #constraints
        //tol1 - convergence tolerane
        template <typename Energy, typename Gradient, typename Hessian,
        typename ConstraintEq, typename JacobianEq, typename Vector, typename Direction, typename Linesearch>
        bool optimizeWithLineSearch(Vector &x0, Energy &f, Gradient &g, Hessian &H, ConstraintEq &ceq, JacobianEq &Aeq, Direction &solver, Linesearch &linesearch, double tol1 = 1e-5, unsigned int maxIter = UINT_MAX) {
            
            unsigned int iter = 0;
            
            while(iter < maxIter) {
                iter++;
                
                //compute search direction
                Vector p = solver(H(x0), g(x0), ceq(x0), Aeq(x0), x0);
                
                //line search and check for convergence
                if (linesearch(x0, p, f,g,H, ceq, Aeq, tol1) == true) {
                    break;
                }
                
            }
            
            return true;
            
        }
        
        
        
    }
    
}
#endif //NEWTON_HPP
