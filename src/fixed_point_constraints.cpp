#include <fixed_point_constraints.h>
#include <algorithm>
void fixed_point_constraints(Eigen::SparseMatrixd &P, unsigned int q_size, const std::vector<unsigned int> indices) {

    unsigned int l = indices.size();
    unsigned int n = q_size / 3;
    //std::cout<<3*l<<std::endl;
    std::vector<int> non_indices, vis(n, 0);
    for(int i = 0; i < l; i++)
    {
        int index = indices[i];
        vis[index] = 1;
    }

    for(int i = 0; i < n; i++)
    {
        if(!vis[i])
        {
            non_indices.push_back(i);
        }
    }

    P.resize(q_size - 3*l, q_size);
    P.setZero();

    for(int i = 0; i < non_indices.size(); i++)
    {
        int index = non_indices[i];
        for(int j=0;j<3;j++)
        {
            P.coeffRef(3*i + j, 3*index + j) = 1;
        }
    }

}