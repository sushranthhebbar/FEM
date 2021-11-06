#include <levelset.h>
#include <math.h>
#include <algorithm>
#include <iostream>
void levelset(Eigen::VectorXd &phi, Eigen::Ref<const Eigen::Vector3d> corner, double cell_width, int grid_length, Eigen::Ref<const Eigen::VectorXi> Ib, Eigen::Ref<const Eigen::VectorXd> q){

    int n = phi.rows();
    for(int i = 0; i < n; i++){
        int z = i / (grid_length * grid_length);
        int tmp = i - z * grid_length * grid_length;
        int y = tmp / grid_length;
        //y = grid_length - y;
        int x = tmp % grid_length;
        Eigen::Vector3d new_point;
        new_point << x, y, z;
        new_point *= cell_width;
        new_point+= corner;
        double min_dist = 1e9;
        for(int j = 0; j < Ib.rows(); j++){
            Eigen::Vector3d dist = q.segment<3>(3 * Ib(j)) - new_point;
            double dist_norm = sqrt(dist.transpose() * dist);
            min_dist = std::min(min_dist, dist_norm);
        }
        phi(i) = min_dist;
    }

}