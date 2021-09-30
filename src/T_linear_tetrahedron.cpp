#include <T_linear_tetrahedron.h>
#include <mass_matrix_linear_tetrahedron.h>
#include <iostream>
void T_linear_tetrahedron(double &T, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::RowVectorXi> element, double density, double volume) {

    Eigen::Matrix1212d M;
    mass_matrix_linear_tetrahedron(M, qdot, element, density, volume);
    Eigen::Vector3d q0, q1, q2, q3;
    Eigen::Vector12d q;
    q0 = qdot.segment<3>(3 * element(0));
    q1 = qdot.segment<3>(3 * element(1));
    q2 = qdot.segment<3>(3 * element(2));
    q3 = qdot.segment<3>(3 * element(3));
    q << q0, q1, q2, q3;
    T = 0.5 * q.transpose() * M * q;
}