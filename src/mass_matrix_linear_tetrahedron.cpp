 
 #include <mass_matrix_linear_tetrahedron.h>

 void mass_matrix_linear_tetrahedron(Eigen::Matrix1212d &M, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::RowVectorXi> element, double density, double volume) {
                
    Eigen::MatrixXd A, B, I;
    A.resize(3, 3);
    B.resize(3, 3);
    I.resize(3, 3);
    I.setIdentity();
    A = density * volume * 6 * I / 60;
    B = density * volume * 6 * I / 120;
    M << A, B, B, B, B, A, B, B, B, B, A, B, B, B, B, A;
 }