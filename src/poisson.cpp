#include<poisson.h>
#include<iostream>

void poisson(Eigen::MatrixXd &potential, Eigen::Ref<const Eigen::VectorXd> theta, double cell_width, double k, int grid_length){

    potential.resize(theta.rows(), 3);
    Eigen::SparseMatrixd K(theta.rows(), theta.rows());
    std::vector<std::vector<Eigen::Triplet<double>>> tripletList;
    tripletList.resize(3);
    Eigen::VectorXd B(theta.rows());
    B.setZero();
    for(int i = 0; i < theta.rows(); i++){
        if(theta(i) == 1){
            for(int j = 0; j < 3; j++){
                tripletList[j].push_back({i, i, 1});
            }
        }
        else{
            int a[] = {0, 0, 0};
            int z = i / (grid_length * grid_length);
            int tmp = i - z * grid_length * grid_length;
            int y = tmp / grid_length;
            //y = grid_length - y;
            int x = tmp % grid_length;
            double sus_i = k * (1 - theta(i));
            for(int j = 0; j < 3; j++){
                a[j] = 1;
                int nxt = grid_length * grid_length * (z + a[2]) + grid_length * (y + a[1]) + x + a[0];
                int prev = grid_length * grid_length * (z - a[2]) + grid_length * (y - a[1]) + x - a[0];
                double sus_nxt = k * (1 - theta(nxt));
                double sus_prev = k * (1 - theta(prev));
                double p = (sus_i + sus_nxt) / 2.0;
                double q = (sus_i + sus_prev) / 2.0;
                tripletList[j].push_back({i , i, -1 * (p + q)});
                if(theta(nxt) != 1){
                    if(nxt>=0 && nxt<=theta.rows()){
                        tripletList[j].push_back({i, nxt, p});
                    }
                }
                if(theta(prev) != 1){
                    if(prev>=0 && prev<=theta.rows()){
                        tripletList[j].push_back({i, prev, q});
                    }
                }
                a[j] = 0;
            }
        }
    }

    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> cg_x, cg_y, cg_z, cg;
    //std::cout<<"Before final for loop in poisson"<<std::endl;
    for(int j = 0; j < 3; j++){
        K.setZero();
        K.setFromTriplets(tripletList[j].begin(), tripletList[j].end());
        cg.compute(K);
        potential.col(j) = cg.solve(B);
    }
    /*K.setZero();
    K.setFromTriplets(tripletList[0].begin(), tripletList[0].end());
    cg_x.compute(K);
    potential.col(0) = cg_x.solve(B);

    K.setZero();
    K.setFromTriplets(tripletList[1].begin(), tripletList[1].end());
    cg_y.compute(K);
    potential.col(1) = cg_y.solve(B);

    K.setZero();
    K.setFromTriplets(tripletList[2].begin(), tripletList[2].end());
    //std::cout<<tripletList[2].size()<<std::endl;
    cg_z.compute(K);
    potential.col(2) = cg_z.solve(B);*/
    
}