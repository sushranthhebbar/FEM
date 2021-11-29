#include<poisson.h>
#include<bar_magnet.h>
#include<iostream>

void poisson(Eigen::MatrixXd &potential, Eigen::Ref<const Eigen::VectorXd> theta, double cell_width, double k, int grid_length, bool constant, Eigen::Ref<const Eigen::Vector3d> corner, double force, Eigen::Ref<const Eigen::MatrixXd> Po1){

    //std::cout<<"HERE"<<std::endl;
    potential.resize(theta.rows(), 3);
    Eigen::SparseMatrixd K(theta.rows(), theta.rows());
    std::vector<std::vector<Eigen::Triplet<double>>> tripletList;
    tripletList.resize(3);
    Eigen::VectorXd C(theta.rows());
    C.setZero();
    std::vector<Eigen::VectorXd> B(3);
    for(int i = 0; i < 3; i++) {
        B[i].resize(theta.rows());
        B[i].setZero();
    }
    for(int i = 0; i < theta.rows(); i++){
        if(theta(i) == 1){
            for(int j = 0; j < 3; j++){
                tripletList[j].push_back({i, i, 1});
            }
        }
        else{
            int a[] = {0, 0, 0};
            double sus_i = k * (1 - theta(i));
            for(int j = 0; j < 3; j++){
                a[j] = 1;
                int z = i / (grid_length * grid_length);
                int tmp = i - z * grid_length * grid_length;
                int y = tmp / grid_length;
                //y = grid_length - y;
                int x = tmp % grid_length;

                int nxt = grid_length * grid_length * (z + a[2]) + grid_length * (y + a[1]) + x + a[0];
                int prev = grid_length * grid_length * (z - a[2]) + grid_length * (y - a[1]) + x - a[0];
                double sus_nxt = k * (1 - theta(nxt));
                double sus_prev = k * (1 - theta(prev));
                double p = (sus_i + sus_nxt) / 2.0;
                double q = (sus_i + sus_prev) / 2.0;
                tripletList[j].push_back({i , i, -1 * (p + q)});
                if(!constant){
                    Eigen::VectorXd H(3);
                    Eigen::Vector3d Po2;
                    Po2 << x, y, z;
                    Po2 *= cell_width;
                    Po2+= corner;
                    bar_magnet(H, Po1, Po2, force);
                    B[j](i) = -1 * (p + q) * H(j) * cell_width;
                }
                if(theta(nxt) != 1){
                    if(nxt>=0 && nxt<=theta.rows()){
                        tripletList[j].push_back({i, nxt, p});
                    }
                    if(!constant){
                        Eigen::VectorXd H(3);
                        Eigen::Vector3d Po2;
                        z = nxt / (grid_length * grid_length);
                        tmp = nxt - z * grid_length * grid_length;
                        y = tmp / grid_length;
                        //y = grid_length - y;
                        x = tmp % grid_length;
                        Po2 << x, y, z;
                        Po2 *= cell_width;
                        Po2+= corner;
                        bar_magnet(H, Po1, Po2, force);
                        B[j](i) += p * H(j) * cell_width;
                    }
                }
                if(theta(prev) != 1){
                    if(prev>=0 && prev<=theta.rows()){
                        tripletList[j].push_back({i, prev, q});
                    }
                    if(!constant){
                        Eigen::VectorXd H(3);
                        Eigen::Vector3d Po2;
                        z = prev / (grid_length * grid_length);
                        tmp = prev - z * grid_length * grid_length;
                        y = tmp / grid_length;
                        //y = grid_length - y;
                        x = tmp % grid_length;
                        Po2 << x, y, z;
                        Po2 *= cell_width;
                        Po2+= corner;
                        bar_magnet(H, Po1, Po2, force);
                        B[j](i) += q * H(j) * cell_width;
                    }
                }
                a[j] = 0;
            }
        }
    }

    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> cg;
    //std::cout<<"Before final for loop in poisson"<<std::endl;
    if(!constant){
            for(int j = 0; j < 3; j++){
                K.setZero();
                K.setFromTriplets(tripletList[j].begin(), tripletList[j].end());
                cg.compute(K);
                potential.col(j) = cg.solve(B[j]);
            }
    }
    else{
            for(int j = 0; j < 3; j++){
                K.setZero();
                K.setFromTriplets(tripletList[j].begin(), tripletList[j].end());
                cg.compute(K);
                potential.col(j) = cg.solve(C);
            }
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