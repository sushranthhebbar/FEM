#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>

void poisson(Eigen::MatrixXd &potential, Eigen::Ref<const Eigen::VectorXd> theta, double cell_width, double k, int grid_length, bool constant, Eigen::Ref<const Eigen::Vector3d> corner, double force, Eigen::Ref<const Eigen::MatrixXd> Po1);