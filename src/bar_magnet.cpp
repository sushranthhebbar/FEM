#include<bar_magnet.h>
#include<iostream>
void bar_magnet(Eigen::VectorXd &H, Eigen::Ref<const Eigen::MatrixXd> Po1, Eigen::Ref<const Eigen::VectorXd> Po2, double force){
    //std::cout<<"Inside bar magent"<<std::endl;
    Eigen::Vector3d r = Po1.row(0) - Po2;
    double f = force / (r.transpose() * r);
    double phi = atan(sqrt(r(0) * r(0) + r(1) * r(1))/ r(2));
    double theta = atan(r(1) / r(0));
    H(0) = f * cos(phi);
    H(1) = f * sin(phi) * cos(theta);
    H(2) = f * sin(phi) * sin(theta);
}