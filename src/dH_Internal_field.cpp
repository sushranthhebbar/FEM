#include <dH_Internal_field.h>
#include <iostream>

void dH_Internal_field(Eigen::MatrixXd &dH, Eigen::Ref<const Eigen::MatrixXd> potential, double cell_width, int grid_length){

    dH.resize(potential.rows(), potential.cols());
    //std::cout<<potential.rows()<<std::endl;
    for(int i = 0; i < potential.rows(); i++){
        int a[] = {0, 0, 0};
        int z = i / (grid_length * grid_length);
        int tmp = i - z * grid_length * grid_length;
        int y = tmp / grid_length;
        //y = grid_length - y;
        int x = tmp % grid_length;
        int b[] = {x, y, z};
        for(int j = 0; j < 3; j++){
            a[j] = 1;
            int nxt = grid_length * grid_length * (z + a[2]) + grid_length * (y + a[1]) + x + a[0];
            int prev = grid_length * grid_length * (z - a[2]) + grid_length * (y - a[1]) + x - a[0];
            double p = potential(nxt, j);
            double q = potential(prev, j);
            if(b[j] == grid_length - 1){
                p = 0.0;
            }
            if(b[j] == 0){
                q = 0.0;
            }
            dH(i, j) = (p - q)/(2 * cell_width);
            a[j] = 0;
        }
    }
}