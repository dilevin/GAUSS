#ifndef ANDERSON_HPP
#define ANDERSON_HPP
#include <vector>
#include <climits>


namespace Gauss {
    namespace Optimization {

        // Comment: dissapointingly this gives me no speed up for newton's solve using a direct solver. More experiments required but my guess is that the effort spent
        //to compute the search direction is to great and the direction to good locally such that there's nothing for anderson acceleration to help with.
        //might also be getting killed by assembly as well.
        
        //Outer loop for a linesearch based optimization
        //  Energy - a function/functor of the form f(x) that returns the energy associated with the  iterate x
        //  ConstraintEq - a function/functor of the form C(x) that returns the vector valued constraint function associated with iterate x
        // JacobianEq -   a function/functor of the form J(x) that returns the gradient of C(x) associated with the iterate x
        //F is an nxm matrix where n is the number of DOFs in the system and is the size of the anderson acceleration basis
        //dF is an nxm-1 matrix for storing the F deltas required for the anderson acceleration step
        // x0 - initial point for the solver. For constrained solves, x0 is the stacked vector of primal and dual variables. Must be of size #primal + #constraints
        //tol1 - convergence tolerane
        template <typename Energy, typename Gradient, typename Hessian,
        typename ConstraintEq, typename JacobianEq, typename Vector, typename Matrix, typename Solver, typename LineSearch, typename Callback>
        bool optimizeAndersonAcceleration(Vector &x0, Energy &f, Gradient &g, Hessian &H, ConstraintEq &ceq, JacobianEq &Aeq, Matrix &F, Matrix &dF, Matrix &G, Matrix &dG, Solver &solver, LineSearch &linesearch, Callback &scallback, double tol1 = 1e-5, unsigned int maxIter = UINT_MAX) {
            
            auto J = *Aeq(x0);
            auto b = *ceq(x0);
            Vector q;
            Vector p;
            Vector theta;
            
            //Vector p = solver(H(x0), gradient, ceq(x0), Aeq(x0), x0);
            unsigned int m = F.cols();
            unsigned int iter = 0;
        
            //compute initial F (direction) and Energy (E0)
            double E0,E1;
            F.col(0) = solver(H(x0), g(x0), ceq(x0), Aeq(x0), x0);
            G.col(0) = x0 + F.col(0);
            q = x0;
            E0 = f(q) + 100000.0*(J*x0.head(J.cols()) - b).norm();
            
            while(iter < maxIter) {
                
                iter++;
                unsigned int k = std::min(m-1, iter);
                
                //make sure Qk decreases energy
                p = F.col(k-1);
                linesearch(q, p, f,g,H, ceq, Aeq, scallback, tol1);
                
                scallback(q);
                J = *Aeq(q);
                b = *ceq(q);
                E1 = f(q) + 100000.0*(J*G.col(k).head(J.cols()) - b).norm();
                
                if(E1 > E0) {
                    p = F.col(k-1);
                    q = G.col(k-1) - F.col(k-1);
                    linesearch(q, p, f,g,H, ceq, Aeq, scallback, tol1);
                    
                    scallback(q);
                    J = *Aeq(q);
                    b = *ceq(q);
                    E1 = f(q) + 100000.0*(J*G.col(k).head(J.cols()) - b).norm();
                }
                
                if(fabs(E0-E1)/fabs(E0) < tol1) {
                    x0 = q;
                    return true;
                }
                
                E0 = E1;
                
                F.col(k) = solver(H(q), g(q), ceq(q), Aeq(q), q);
                G.col(k) = q + F.col(k);
                
                for(unsigned int kk=0; kk < k; ++kk) {
                    dF.col(kk) = F.col(kk) - F.col(k);
                    dG.col(kk) = G.col(kk) - G.col(k);
                }
                
                //compute thetas (alg 1, line 10)
                theta = (dF.leftCols(k).transpose()*dF.leftCols(k)).colPivHouseholderQr().solve(dF.leftCols(k).transpose()*F.col(k));
                
                q = G.col(k) - dG.leftCols(k)*theta;
                
                scallback(q);
                //E1 = f(q) + 100000.0*(J*q.head(J.cols()) - b).norm();
                
                //make storage space
                if((iter+1) >= m) {
                    F.leftCols(m-1) = F.rightCols(m-1);
                    G.leftCols(m-1) = G.rightCols(m-1);
                }
            }
            
            return true;
            
        }
    }
}

#endif ANDERSON_HPP
