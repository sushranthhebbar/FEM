#include <interpolate.h>
#include <iostream>

void interpolate(Eigen::VectorXd &boundary_value, Eigen::Ref<const Eigen::MatrixXd> dH, Eigen::Ref<const Eigen::Vector3d> corner, double cell_width, Eigen::Ref<const Eigen::VectorXi> Ib, Eigen::Ref<const Eigen::VectorXd> q, int grid_length){

    boundary_value.resize(3  * Ib.rows());
    boundary_value.setZero();

    for(int i = 0; i < Ib.rows(); i++){
        Eigen::Vector3d point= q.segment<3>(3 * Ib(i));
        Eigen::Vector3d diff = (point - corner) / cell_width;
        int i000 = floor(diff(2)) * grid_length * grid_length + floor(diff(1)) * grid_length + floor(diff(0));
        int i100 = i000 + 1;
        int i010 = i000 + grid_length;
        int i110 = i000 + 1 + grid_length;
        int i001 = i000 + grid_length * grid_length;
        int i011 = i000 + grid_length * grid_length + grid_length;
        int i101 = i000 + 1 + grid_length * grid_length;
        int i111 = i000 + 1 + grid_length * grid_length + grid_length;

        Eigen::Vector3d bl;
        bl << floor(diff(0)), floor(diff(1)), floor(diff(2));
        bl = bl * cell_width;
        bl = corner + bl;

        Eigen::Vector3d coeff = (point - bl) / cell_width;

        Eigen::Vector3d c00 = (1 - coeff(0)) * dH.row(i000) + coeff(0) * dH.row(i100);
        Eigen::Vector3d c01 = (1 - coeff(0)) * dH.row(i001) + coeff(0) * dH.row(i101);
        Eigen::Vector3d c10 = (1 - coeff(0)) * dH.row(i010) + coeff(0) * dH.row(i110);
        Eigen::Vector3d c11 = (1 - coeff(0)) * dH.row(i011) + coeff(0) * dH.row(i111);

        Eigen::Vector3d c0 = (1 - coeff(1)) * c00 + coeff(1) * c10;
        Eigen::Vector3d c1 = (1 - coeff(1)) * c01 + coeff(1) * c11;

        Eigen::Vector3d c = (1 - coeff(2)) * c0 + coeff(2) * c1;

        boundary_value.segment<3>(3 * i)  = c;

    }
}