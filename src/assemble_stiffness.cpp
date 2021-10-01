#include <assemble_stiffness.h>

void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T, Eigen::Ref<const Eigen::VectorXd> v0, 
                     double C, double D) { 
        
        
        std::vector<Eigen::Triplet<double>> tripletList;
        Eigen::Matrix1212d H;
        for(int i = 0; i < T.rows(); i++)
        {
            Eigen::RowVectorXi element = T.row(i);
            d2V_linear_tetrahedron_dq2(H, q, V, element, v0(i), C, D);

            for(int j = 0; j < 4; j++)
            {
                for(int k = 0; k < 4; k++)
                {
                    for(int ii = 0; ii < 3; ii++)
                    {
                        for(int jj = 0; jj < 3; jj++)
                        {
                            Eigen::Triplet<double> t  = {3 * element(j) + ii, 3 * element(k) + jj, -1 * H(3 * j + ii, 3 * k + jj)};
                            tripletList.push_back(t);
                        }
                    }
                }
            }
        }

        K.resize(q.size(), q.size());
        K.setFromTriplets(tripletList.begin(), tripletList.end());

    };
