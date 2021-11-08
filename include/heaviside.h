#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>

void heaviside(Eigen::VectorXd &theta, Eigen::Ref<const Eigen::VectorXd> phi, double epsilon);