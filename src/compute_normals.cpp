#include <compute_normals.h>
#include <iostream>
void compute_normals(Eigen::MatrixXd &N ,Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXi> Ib, Eigen::Ref<const Eigen::MatrixXi> Fb){

    //std::cout<<"Inside compute_normals"<<std::endl;
    //N.resize(Fb.rows(), Fb.cols());
    Eigen::MatrixXd Vb(Ib.rows(), Fb.cols());
    //std::cout<<"Before for loop"<<std::endl;
    for(int i = 0; i < Ib.rows(); i++)
    {
        int idx = Ib(i);
        //std::cout<<i<<std::endl;
        Eigen::Vector3d r = q.segment<3>(3 * idx);
        //std::cout<<r<<std::endl;
        Vb.row(i) = r.transpose();
    }
    //std::cout<<"Before per vertex"<<std::endl;
    //std::cout<<Vb<<std::endl;
    //std::cout<<std::endl;
    //igl::per_vertex_normals(Vb, Fb, igl::PER_VERTEX_NORMALS_WEIGHTING_TYPE_UNIFORM, N);
    igl::per_vertex_normals(Vb, Fb, N);
    //std::cout<<N.rows()<<std::endl;

}