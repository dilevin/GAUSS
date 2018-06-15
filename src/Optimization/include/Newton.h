#ifndef NEWTON_HPP
#define NEWTON_HPP
#include <vector>

namespace Gauss {
    namespace Optimization {
        
        //Test Functions, Delete in a bit
        void printVal(double val) { std::cout<<val<<"\n"; }
        
        template<typename A>
        void doSomething(A a) {
            std::cout<<a<<"\n";
        }
    
        //functions to implement an equality constrained newton solve using Gauss
        //Energy
        template<typename Vector, typename World, typename MassMatrix, typename QDot>
        inline auto & energyImplicitWorld(Vector &x0, World &world, MassMatrix &massMatrix, QDot &qDot) {
            return (getEnergy(world) - mapStateEigen<1>(world).transpose()*massMatrix*qDot);
        }
        
        //Gradient
        template<typename Vector, typename World, typename MassMatrix, typename ForceVector, typename QDot, typename DataType>
        inline auto & gradientImplicit(Vector &x0, World &world, MassMatrix &massMatrix, ForceVector &forceVector, DataType &dt, QDot &qDot) {
            ASSEMBLEVECINIT(forceVector, world.getNumQDotDOFs());
            ASSEMBLELIST(forceVector, world.getForceList(), getForce);
            ASSEMBLELIST(forceVector, world.getSystemList(), getForce);
            ASSEMBLEEND(forceVector);
            
            (*forceVector) *= -dt;
            (*forceVector) += massMatrix*(mapStateEigen<1>(world)-qDot);
            return (*forceVector);
        };
        
        //Hessian actually builds the KKT for world stuff (it's faster than merging sparse matrices)
        template<typename Vector, typename World, typename MassMatrix, typename StiffnessMatrix, typename ConstraintEq, typename QDot, typename DataType>
        inline auto & hessianImplicit(Vector &x0, World &world,
                                      MassMatrix &massMatrix, StiffnessMatrix &stiffnessMatrix, ConstraintEq &ceq,
                                      DataType &dt) {
            //get stiffness matrix
            ASSEMBLEMATINIT(stiffnessMatrix, world.getNumQDotDOFs()+world.getNumConstraints(), world.getNumQDotDOFs()+world.getNumConstraints());
            ASSEMBLELIST(stiffnessMatrix, world.getSystemList(), getStiffnessMatrix);
            ASSEMBLELIST(stiffnessMatrix, world.getForceList(), getStiffnessMatrix);
            
            ASSEMBLEEND(stiffnessMatrix);
            
            (*stiffnessMatrix) *= -(dt*dt);
            (*stiffnessMatrix) += massMatrix;
            return (*stiffnessMatrix);
        };
                                        
                                        
        //Equality Constraints (velocity level)
        //Hessian actually builds the KKT for world stuff (it's faster than merging sparse matrices)
        template<typename Vector, typename World, typename MassMatrix, typename StiffnessMatrix, typename ConstraintEq, typename QDot, typename DataType>
        inline auto & hessianImplicit(Vector &x0, World &world,
                                      MassMatrix &massMatrix, StiffnessMatrix &stiffnessMatrix, ConstraintEq &ceq,
                                      DataType &dt) {
            //get stiffness matrix
            ASSEMBLEMATINIT(stiffnessMatrix, world.getNumQDotDOFs()+world.getNumConstraints(), world.getNumQDotDOFs()+world.getNumConstraints());
            ASSEMBLELIST(stiffnessMatrix, world.getSystemList(), getStiffnessMatrix);
            ASSEMBLELIST(stiffnessMatrix, world.getForceList(), getStiffnessMatrix);
            
            ASSEMBLEEND(stiffnessMatrix);
            
            (*stiffnessMatrix) *= -(dt*dt);
            (*stiffnessMatrix) += massMatrix;
            return (*stiffnessMatrix);
        };
        
                                        
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
        inline bool backTracking(Energy &f, ConstraintEq &ceq, JacobianEq &Aeq, Direction &solver, Vector &x0, StepCallback &scallback, double tol1 = 1e-5) {
            
            //std::cout<<"NORM2: "<<x0.norm()<<"\n";
            //first version is a lazy Newton step with no line search.
            
            scallback(x0);
            double E0 = f(x0);
            Vector p = solver(H(x0), g(x0), Aeq(x0), ceq(x0));
            double gStep = g(x0).transpose()*p;
            
            //back tracking line search
            double alpha =1;
            double c = 1e-4; //from Nocedal and Wright pg 31
            double rho = 0.5;
            
            while(true) {
                
                scallback(x0+alpha*p);
                
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
        
        //Outer loop for a linesearch based optimization
        //  Energy - a function/functor of the form f(x) that returns the energy associated with the  iterate x
        //  ConstraintEq - a function/functor of the form C(x) that returns the vector valued constraint function associated with iterate x
        // JacobianEq -   a function/functor of the form J(x) that returns the gradient of C(x) associated with the iterate x
        // LineSearch - A function/function of the form  linesearch(f,g,H, ceq, Aeq, x0, tol1) that computes a search direction then does linesearch
        // x0 - initial point for the solver. For constrained solves, x0 is the stacked vector of primal and dual variables. Must be of size #primal + #constraints
        //tol1 - convergence tolerane
        template <typename Energy, typename Gradient, typename Hessian,
        typename ConstraintEq, typename JacobianEq,
        typename StepCallback, typename Vector, typename Linesearch>
        bool optimizeWithLineSearch(Energy &f, ConstraintEq &ceq, JacobianEq &Aeq, LineSearch &linesearch, Vector &x0, double tol1 = 1e-5) {
            
            unsigned int iter = 0;
            
            while(iter < maxIter && linesearch(f,g,H, ceq, Aeq, x0, tol1) == false) {
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
