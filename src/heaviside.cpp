#include <heaviside.h>
#include <iostream>
#include <algorithm>

void heaviside(Eigen::VectorXd &theta, Eigen::Ref<const Eigen::VectorXd> phi, double epsilon){

    double pi = 3.1415;
    theta.resize(phi.size());
    for(int i = 0; i < phi.rows(); i++){
        double phi_i = phi(i);
        if(phi_i <= -1 * epsilon){
            theta(i) = 0;
        }
        else if(phi_i >= epsilon){
            theta(i) = 1;
        }
        else{
            theta(i) = 0.5 + 0.5 * phi_i / epsilon + 0.5 * sin(pi * phi_i / epsilon) / pi;
        }
    }

}