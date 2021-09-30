#include <dV_linear_tetrahedron_dq.h>

#include <dphi_linear_tetrahedron_dX.h>
#include <dpsi_neo_hookean_dF.h>
#include <quadrature_single_point.h>
#include <iostream>

void dV_linear_tetrahedron_dq(Eigen::Vector12d &dV, Eigen::Ref<const Eigen::VectorXd> q, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                          double C, double D) {

   auto neohookean_linear_tet = [&](Eigen::Vector12d &dV, Eigen::Ref<const Eigen::VectorXd>q, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
        
       Eigen::Matrix34d q4; 
       Eigen::Matrix43d dphi; 
       Eigen::Vector3d uX;
       Eigen::Matrix3Xd F;
       Eigen::Vector9d dF;
       Eigen::MatrixXd B;
       B.resize(12, 9);
       for(int i = 0; i < element.size(); i++)
       {
           Eigen::Vector3d qi = q.segment<3>(3 * element(i));
           q4.col(i) = qi;
           uX = uX + qi * ((i != 3) ? X(i) : 1 - X(0) - X(1) - X(2));
       }

       dphi_linear_tetrahedron_dX(dphi, V, element, uX);
       F = q4 * dphi;
       dpsi_neo_hookean_dF(dF, F, C, D);

       Eigen::MatrixXd A1, A2, A3;
       A1.resize(12, 4);
       A2.resize(12, 4);
       A3.resize(12, 4);
       A1.setZero();
       A2.setZero();
       A3.setZero();
       for(int i = 0; i < 3; i++)
       {
           for(int j = 0; j < 4; j++)
           {
               if(i==0)
               {
                    A1(3*j + i,j) = 1;
               }
               else if(i==1)
               {
                    A2(3*j + i,j) = 1;
               }
               else
               {
                    A3(3*j + i,j) = 1;
               }
           }
       }
       B << (A1 * dphi), (A2 * dphi), (A3 * dphi);
       dV = B * dF;

    };

    quadrature_single_point(dV, q, element, volume, neohookean_linear_tet);  
    
}