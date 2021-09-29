#include <Eigen/Dense>
#include <EigenTypes.h>
#include <iostream>
using namespace std;

int main()
{
    Eigen::VectorXd vec1(3);
    vec1 << 1, 0, 0;

    Eigen::VectorXd vec2(3);
    vec2 << 4, 5, 6;

    Eigen::VectorXd vec3(3);
    vec3 << 7, 8, 9;

    //cout<<"HERE"<<endl;

    Eigen::MatrixXd m(3, 3);
    m << vec1, vec2, vec3;
    Eigen::Matrix34d phi;
    /*phi(0) = 1;
    phi(1) = vec3(0);
    phi(2) = vec3(1);
    phi(3) = vec3(2);*/
    //phi << (1), vec3;

    //phi << (1);
    //phi << vec3;
    //Eigen::Ref<const Eigen::Vector3d> x = {1 , 2 , 3};
    //cout<<phi<<endl;
    //cout<<m<<endl;
    //cout<<(m.inverse() * vec1)<<endl;

    Eigen::Matrix43d dphi;
    Eigen::Vector3d one;
    one << 1, 1, 1;
    //cout<<one.transpose()<<endl;
    //dphi << (one.transpose()) , m;
    //dphi << m;
    //cout<<dphi<<endl;

    /*phi.col(0) = vec1;
    phi.col(1) = vec2;
    phi.col(2) = vec3;
    phi.col(3) = one;
    cout<<phi<<endl;*/

    Eigen::Vector6d q;
    q << 1, 2, 3, 4, 5, 6;
    Eigen::Vector3d qi = q.segment<3>(2);
    cout<<qi<<endl;

    return 0;
}