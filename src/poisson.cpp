#include<poisson.h>
#include<iostream>

void poisson(Eigen::VectorXd &potential, Eigen::Ref<const Eigen::VectorXd> theta, double cell_width, double k, int grid_length){

    potential.resize(theta.rows());
    Eigen::SparseMatrixd K(theta.rows(), theta.rows());
    std::vector<Eigen::Triplet<double>> tripletList;
    Eigen::VectorXd B(theta.rows());
    B.setZero();
    for(int i = 0; i < theta.rows(); i++){
        if(theta(i) == 1){
            tripletList.push_back({i, i, 1});
        }
        else{
            int z = i / (grid_length * grid_length);
            int tmp = i - z * grid_length * grid_length;
            int y = tmp / grid_length;
            //y = grid_length - y;
            int x = tmp % grid_length;
            int nxt = grid_length * grid_length * z + grid_length * y + x + 1;
            int prev = grid_length * grid_length * z + grid_length * y + x - 1;
            double sus_nxt = k * (1 - theta(nxt));
            double sus_prev = k * (1 - theta(prev));
            double sus_i = k * (1 - theta(i));
            double p = (sus_i + sus_nxt) / 2.0;
            double q = (sus_i + sus_prev) / 2.0;
            tripletList.push_back({i , i, -1 * (p + q)});
            if(theta(nxt) != 1){
                tripletList.push_back({i, nxt, p});
            }
            if(theta(prev) != 1){
                tripletList.push_back({i, prev, q});
            }
        }
    }
    K.setFromTriplets(tripletList.begin(), tripletList.end());
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> cg;
    cg.compute(K);
    potential = cg.solve(B);
}