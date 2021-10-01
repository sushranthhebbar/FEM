#include <Eigen/Dense>
#include <Eigen/Sparse>
#include<Eigen/SparseCholesky>
#include <EigenTypes.h>

//Input:
//  q - generalized coordinates for the FEM system
//  qdot - generalized velocity for the FEM system
//  dt - the time step in seconds
//  mass - the mass matrix
//  force(f, q, qdot) - a function that computes the force acting on the FEM system. This takes q and qdot as parameters, returns the force in f.
//  stiffness(K, q, qdot) - a function that computes the stiffness (negative second derivative of the potential energy). This takes q and qdot as parameters, returns the stiffness matrix in K.  
//  tmp_force - scratch space to collect forces
//  tmp_stiffness - scratch space to collect stiffness matrix
//Output:
//  q - set q to the updated generalized coordinate using linearly implicit time integration
//  qdot - set qdot to the updated generalized velocity using linearly implicit time integration
template<typename FORCE, typename STIFFNESS> 
inline void linearly_implicit_euler(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, 
                            const Eigen::SparseMatrixd &mass,  FORCE &force, STIFFNESS &stiffness, 
                            Eigen::VectorXd &tmp_force, Eigen::SparseMatrixd &tmp_stiffness) {
    
    force(tmp_force, q, qdot);
    stiffness(tmp_stiffness, q, qdot);

    tmp_force = mass * qdot + dt * tmp_force;
    tmp_stiffness = mass - dt * dt * tmp_stiffness;

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(tmp_stiffness);
    if(solver.info()!=Eigen::Success) 
    {
        std::cout<<solver.info()<<std::endl;
        // decomposition failed
        return;
    }
    qdot = solver.solve(tmp_force);
    if(solver.info()!=Eigen::Success)
    {
        std::cout<<solver.info()<<std::endl;
        // solving failed
        return;
    }

    q = q + qdot * dt;
}