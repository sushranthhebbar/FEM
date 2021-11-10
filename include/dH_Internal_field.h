#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>

void dH_Internal_field(Eigen::MatrixXd &dH, Eigen::Ref<const Eigen::MatrixXd> potential, double cell_width, int grid_length);