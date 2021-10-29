#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>

void compute_normals(Eigen::MatrixXd &N ,Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXi> Vb, Eigen::Ref<const Eigen::MatrixXi> Fb);