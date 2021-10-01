#include <mass_matrix_mesh.h>
#include <mass_matrix_linear_tetrahedron.h>

void mass_matrix_mesh(Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::MatrixXi> T, double density, Eigen::Ref<const Eigen::VectorXd> v0) {

        std::vector<Eigen::Triplet<double>> tripletList;
        Eigen::Matrix1212d Mi;
        for(int i = 0; i < T.rows(); i++)
        {
            Eigen::RowVectorXi element = T.row(i);
            mass_matrix_linear_tetrahedron(Mi, qdot, element, density, v0(i));

            for(int j = 0; j < 4; j++)
            {
                for(int k = 0; k < 4; k++)
                {
                    for(int ii = 0; ii < 3; ii++)
                    {
                        for(int jj = 0; jj < 3; jj++)
                        {
                            Eigen::Triplet<double> t  = {3 * element(j) + ii, 3 * element(k) + jj, Mi(3 * j + ii, 3 * k + jj)};
                            tripletList.push_back(t);
                        }
                    }
                }
            }
        }

        M.resize(qdot.size(), qdot.size());
        M.setFromTriplets(tripletList.begin(), tripletList.end());

}