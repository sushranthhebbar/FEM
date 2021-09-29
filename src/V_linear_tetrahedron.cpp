#include <V_linear_tetrahedron.h>

#include <dphi_linear_tetrahedron_dX.h>
#include <psi_neo_hookean.h>
#include <quadrature_single_point.h>

void V_linear_tetrahedron(double &energy, Eigen::Ref<const Eigen::VectorXd> q, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                          double C, double D) {

    auto neohookean_linear_tet = [&](double &e, Eigen::Ref<const Eigen::VectorXd>q, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
        
       Eigen::Matrix34d q4; 
       Eigen::Matrix43d dphi; 
       Eigen::Vector3d uX;
       Eigen::Matrix3Xd F;
       for(int i = 0; i < element.size(); i++)
       {
           Eigen::Vector3d qi = q.segment<3>(3 * element(i));
           q4.col(i) = qi;
           uX = uX + qi * ((i != 3) ? X(i) : 1 - X(0) - X(1) - X(2));
       }
       
       dphi_linear_tetrahedron_dX(dphi, V, element, uX);
       F = q4 * dphi;
       psi_neo_hookean(e, F, C, D);
    };

    quadrature_single_point(energy, q, element, volume, neohookean_linear_tet);  
    
}