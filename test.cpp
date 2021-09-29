#include <Eigen/Dense>
#include <EigenTypes.h>
#include <iostream>
using namespace std;

int main()
{
    Eigen::RowVectorXd vec1(3);
    vec1 << 1, 2, 3;

    Eigen::RowVectorXd vec2(3);
    vec2 << 4, 5, 6;

    Eigen::RowVectorXd vec3(3);
    vec3 << 7, 8, 9;


    Eigen::MatrixXd m(3, 3);
    m << vec1, vec2, vec3;

    cout<<m.transpose()<<endl;

    return 0;
}