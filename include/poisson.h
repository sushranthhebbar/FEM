#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>

void poisson(Eigen::VectorXd &potential, Eigen::Ref<const Eigen::VectorXd> theta, double cell_width, double k, int grid_length);