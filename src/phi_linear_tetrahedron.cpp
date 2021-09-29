#include <phi_linear_tetrahedron.h>

void phi_linear_tetrahedron(Eigen::Vector4d &phi, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> x) {

    //Eigen::Ref<const Eigen::Vector3d> X0 = V.row(element(0));
    Eigen::Vector3d X0 = V.row(element(0));
    Eigen::Vector3d X1 = V.row(element(1));
    Eigen::Vector3d X2 = V.row(element(2));
    Eigen::Vector3d X3 = V.row(element(3));

    X1 = X1 - X0;
    X2 = X2 - X0;
    X3 = X3 - X0;
    X0 = x - X0;
    Eigen::MatrixXd T(3, 3);

    T << X1, X2, X3;
    Eigen::Vector3d tmp = T.inverse() * X0;

    phi << (1 - tmp(0) - tmp(1) - tmp(2)), tmp;
}