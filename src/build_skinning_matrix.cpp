#include <build_skinning_matrix.h>
#include <phi_linear_tetrahedron.h>
#include <vector>
#include <iostream>
#include <math.h>

bool SameSide(Eigen::Vector3d v1, Eigen::Vector3d v2, Eigen::Vector3d v3, Eigen::Vector3d v4, Eigen::Vector3d p)
{
    Eigen::Vector3d normal = (v2 - v1).cross(v3 - v1);
    double dotV4 = normal.dot(v4 - v1);
    double dotP = normal.dot(p - v1);
    return signbit(dotV4) == signbit(dotP);
}

bool PointInTetrahedron(Eigen::Vector3d v1, Eigen::Vector3d v2, Eigen::Vector3d v3, Eigen::Vector3d v4, Eigen::Vector3d p)
{
    return SameSide(v1, v2, v3, v4, p) &&
           SameSide(v2, v3, v4, v1, p) &&
           SameSide(v3, v4, v1, v2, p) &&
           SameSide(v4, v1, v2, v3, p);               
}

void build_skinning_matrix(Eigen::SparseMatrixd &N, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T, 
                                                   Eigen::Ref<const Eigen::MatrixXd> V_skin) {

            std::vector<Eigen::Triplet<double>> tripletList;
            for(int i = 0; i < V_skin.rows(); i++)
            {
                Eigen::Vector4d phi;
                int id = -1;
                Eigen::Vector3d x = V_skin.row(i).transpose();
                for(int j = 0; j < T.rows(); j++)
                {
                    Eigen::RowVectorXi element = T.row(j);
                    Eigen::Vector3d v1 = V.row(element(0)).transpose();
                    Eigen::Vector3d v2 = V.row(element(1)).transpose();
                    Eigen::Vector3d v3 = V.row(element(2)).transpose();
                    Eigen::Vector3d v4 = V.row(element(3)).transpose();
                    if(PointInTetrahedron(v1, v2, v3, v4, x)) 
                    {
                        id = j;
                        break;
                    }
                }
                Eigen::RowVectorXi element = T.row(id);
                phi_linear_tetrahedron(phi, V, element, x);
                for(int j = 0; j < 4; j++)
                {
                    Eigen::Triplet<double> t = {i, element(j), phi(j)};
                    tripletList.push_back(t);
                }
            }
            N.resize(V_skin.rows(), V.rows());
            N.setFromTriplets(tripletList.begin(), tripletList.end());
}