#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>
#include <igl/per_vertex_normals.h>

void compute_normals(Eigen::MatrixXd &N ,Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXi> Ib, Eigen::Ref<const Eigen::MatrixXi> Fb);