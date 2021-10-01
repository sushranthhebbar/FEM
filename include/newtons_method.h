#include <Eigen/Dense>
#include <Eigen/Sparse>
#include<Eigen/SparseCholesky>
#include <EigenTypes.h>

//Input:
//  x0 - initial point for newtons search
//  f(x) - function that evaluates and returns the cost function at x
//  g(dx, x) - function that evaluates and returns the gradient of the cost function in dx
//  H(dH, x) - function that evaluates and returns the Hessian in dH (as a sparse matrix).
//  max steps - the maximum newton iterations to take
//  tmp_g and tmp_H are scratch space to store gradients and hessians
//Output: 
//  x0 - update x0 to new value
template<typename Objective, typename Jacobian, typename Hessian>
double newtons_method(Eigen::VectorXd &x0, Objective &f, Jacobian &g, Hessian &H, unsigned int maxSteps, Eigen::VectorXd &tmp_g, Eigen::SparseMatrixd &tmp_H) {
   
   double c = 1e-8;
   for(int i = 0; i < maxSteps; i++)
   {
      g(tmp_g, x0);
      if(tmp_g.norm() < (1e-5)) break;
      H(tmp_H, x0);
      Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
      solver.compute(tmp_H);
      if(solver.info()!=Eigen::Success) 
      {
         std::cout<<solver.info()<<std::endl;
         // decomposition failed
         return 0.0;
      }
      Eigen::VectorXd d = solver.solve(-1 * tmp_g);
      if(solver.info()!=Eigen::Success)
      {
         std::cout<<solver.info()<<std::endl;
         // solving failed
         return 0.0;
      }
      double alpha = 1;
      double lim = f(x0) + c * d.transpose() * tmp_g;
      while(true)
      {
         if((alpha < 1e-5) || (f(x0 + alpha * d) <= lim))
         {
            break;
         }
         alpha = 0.5 * alpha;
      }
      x0 = x0 + alpha * d;
   }

   return 0.0;
}
