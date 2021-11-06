#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>

void levelset(Eigen::VectorXd &phi, Eigen::Ref<const Eigen::Vector3d> corner, double cell_width, int grid_length, Eigen::Ref<const Eigen::VectorXi> Ib, Eigen::Ref<const Eigen::VectorXd> q);