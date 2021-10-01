#include <dV_spring_particle_particle_dq.h>

void dV_spring_particle_particle_dq(Eigen::Ref<Eigen::Vector6d> f, Eigen::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d>     q1, double l0, double stiffness) {

    Eigen::Vector6d q;
    Eigen::Vector3d dq;
    q << q0, q1;
    Eigen:: MatrixXd B(6,6);
    B.setZero();
    for(int i = 0;i < 3; i++)
    {
        B(i, i) = 1;
        B(i+3, i) = -1;
        B(i+3, i+3) = 1;
        B(i, i+3) = -1;
    }
    dq = q1 - q0;
    double c = sqrt(dq.transpose() * dq); 
    c = stiffness * (c - l0) / c;
    f = B * q;
    f = c * f;
    
}