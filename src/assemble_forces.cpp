#include <assemble_forces.h>
#include <iostream>

void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T, Eigen::Ref<const Eigen::VectorXd> v0,
                     double C, double D) { 
        
        f.resize(3 * V.rows());
        f.setZero();
        for(int i = 0; i < T.rows(); i++)
        {
            Eigen::RowVectorXi element = T.row(i);
            Eigen::Vector12d fi;
            dV_linear_tetrahedron_dq(fi, q, V, element, v0(i), C, D);
            fi = -1 * fi;
            for(int i = 0; i < 4; i++)
            {
                int index = element(i);
                for(int j = 0; j < 3; j++)
                {
                    f(3*index + j) = f(3*index + j) + fi(3*i + j);
                }
            }
        }

    };