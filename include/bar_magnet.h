#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>

void bar_magnet(Eigen::VectorXd &H, Eigen::Ref<const Eigen::MatrixXd> Po1, Eigen::Ref<const Eigen::VectorXd> Po2, double force);