#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>

void interpolate(Eigen::VectorXd &boundary_value, Eigen::Ref<const Eigen::MatrixXd> dH, Eigen::Ref<const Eigen::Vector3d> corner, double cell_width, Eigen::Ref<const Eigen::VectorXi> Ib, Eigen::Ref<const Eigen::VectorXd> q);