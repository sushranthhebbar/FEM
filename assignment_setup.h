#ifndef ASSIGNMENT_SETUP_H
#define ASSIGNMENT_SETUP_H

#include <igl/readMESH.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/readOFF.h>
#include <read_tetgen.h>
#include <igl/boundary_facets.h>
#include <igl/volume.h>

//assignment files for implementing simulation and interaction
#include <visualization.h>
#include <init_state.h>
#include <find_min_vertices.h>
#include <fixed_point_constraints.h>
#include <dphi_linear_tetrahedron_dX.h>
#include <psi_neo_hookean.h>
#include <quadrature_single_point.h>
#include <T_linear_tetrahedron.h>
#include <V_linear_tetrahedron.h>
#include <V_spring_particle_particle.h>
#include <dV_linear_tetrahedron_dq.h>
#include <dV_spring_particle_particle_dq.h>
#include <d2V_linear_tetrahedron_dq2.h>
#include <mass_matrix_mesh.h>
#include <assemble_forces.h>
#include <assemble_stiffness.h>
#include <linearly_implicit_euler.h>
#include <implicit_euler.h>
#include <build_skinning_matrix.h>
#include <compute_normals.h>
#include <levelset.h>
#include <heaviside.h>
#include <poisson.h>
#include <dH_Internal_field.h>
#include <interpolate.h>
#include <bar_magnet.h>

//Variable for geometry
Eigen::MatrixXd V; //vertices of simulation mesh 
Eigen::MatrixXi T; //faces of simulation mesh
Eigen::MatrixXi F; //faces of simulation mesh
Eigen::MatrixXi Fb;

Eigen::VectorXi Ib;

//variables for skinning
Eigen::MatrixXd V_skin;
Eigen::MatrixXi F_skin;
Eigen::SparseMatrixd N;

Eigen::MatrixXi Pf;
Eigen::MatrixXd Pv_skin;
Eigen::MatrixXi Pf_skin;
Eigen::SparseMatrixd Pn;

//material parameters
double density = 0.1;
double YM = 6e5; //young's modulus
double mu = 0.4; //poissons ratio
double D = 0.5*(YM*mu)/((1.0+mu)*(1.0-2.0*mu));
double C = 0.5*YM/(2.0*(1.0+mu));
double pi = 3.1415;
double mew = 4 * pi * 1e-7;
double k = 0.33; 
double cell_width = 0.0;
double grid_length = 32;
double epsilon = 0.0;
//BC
std::vector<unsigned int> fixed_point_indices;
Eigen::SparseMatrixd P;
Eigen::VectorXd x0; 
Eigen::Vector3d mi;
Eigen::Vector3d mx;
//mass matrix
Eigen::SparseMatrixd M;
Eigen::VectorXd v0;

//scratch memory for assembly
Eigen::VectorXd tmp_qdot;
Eigen::VectorXd tmp_force;
Eigen::SparseMatrixd tmp_stiffness;
Eigen::MatrixXd Po(1,3);
Eigen::MatrixXd Ini(1,3);
std::vector<std::pair<Eigen::Vector3d, unsigned int>> spring_points;

bool visualization = false;
bool skinning_on = true;
bool fully_implicit = false;
bool bunny = false; 
bool magnet = false;
bool box = false;
bool cube86 = false;
bool constant = true;
bool menu = true;
//selection spring
double k_selected = 1e5;
double mov = 0.0;

inline void simulate(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double t) {  

    double V_ele, T_ele, KE,PE;

    spring_points.clear();
    if(menu){
        Visualize::update_parameters(magnet, constant, Ini, Po);
    }
    //std::cout<<magnet<<std::endl;
    //Interaction spring
    Eigen::Vector3d mouse;
    Eigen::Vector6d dV_mouse;
    Eigen::MatrixXd Nor;
    Eigen::VectorXd boundary_value, H(3);
    double k_selected_now = (Visualize::is_mouse_dragging() ? k_selected : 0.);
    double c = 0.0;

    for(unsigned int pickedi = 0; pickedi < Visualize::picked_vertices().size(); pickedi++) {   
        spring_points.push_back(std::make_pair((P.transpose()*q+x0).segment<3>(3*Visualize::picked_vertices()[pickedi]) + Visualize::mouse_drag_world() + Eigen::Vector3d::Constant(1e-6),3*Visualize::picked_vertices()[pickedi]));
    }

    auto energy = [&](Eigen::Ref<const Eigen::VectorXd> qdot_1)->double {
        double E = 0;
        Eigen::VectorXd newq = P.transpose()*(q+dt*qdot_1)+x0;

        for(unsigned int ei=0; ei<T.rows(); ++ei) {
            
            V_linear_tetrahedron(V_ele,newq , V, T.row(ei), v0(ei), C, D);
            E += V_ele;
        }

        for(unsigned int pickedi = 0; pickedi < spring_points.size(); pickedi++) {   
            V_spring_particle_particle(V_ele, spring_points[pickedi].first, newq.segment<3>(spring_points[pickedi].second), 0.0, k_selected_now);
            E += V_ele;
        }

        E += 0.5*(qdot_1 - qdot).transpose()*M*(qdot_1 - qdot);

        return E;
    };

    auto force = [&](Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q2, Eigen::Ref<const Eigen::VectorXd> qdot2) { 
        
            //std::cout<<q.segment<3>(498*3)<<std::endl;
            assemble_forces(f, P.transpose()*q2+x0, P.transpose()*qdot2, V, T, v0, C,D);
            // between 10^5 and 10^7(too high blows up)
            for(unsigned int pickedi = 0; pickedi < spring_points.size(); pickedi++) {
                dV_spring_particle_particle_dq(dV_mouse, spring_points[pickedi].first, (P.transpose()*q2+x0).segment<3>(spring_points[pickedi].second), 0.0, k_selected_now);
                f.segment<3>(3*Visualize::picked_vertices()[pickedi]) -= dV_mouse.segment<3>(3);
            }
            if(magnet){
                //std::cout<<"Before loop"<<std::endl;
                Eigen::Vector3d corner = V.colwise().minCoeff();
                Eigen::VectorXd phi(32 * 32 * 32);
                Eigen::VectorXd theta;
                Eigen::MatrixXd potential, dH;
                //std::cout<<q.segment<3>(498*3)<<std::endl;
                levelset(phi, corner, cell_width, grid_length, Ib, P.transpose()*q2+x0);
                //std::cout<<q.segment<3>(498*3)<<std::endl;
                epsilon = 3 * cell_width / 2;
                heaviside(theta, phi, epsilon);
                poisson(potential, theta, cell_width, k, grid_length);
                dH_Internal_field(dH, potential, cell_width, grid_length);
                //std::cout<<q.segment<3>(498*3)<<std::endl;
                interpolate(boundary_value, dH, corner, cell_width, Ib, P.transpose()*q2+x0, grid_length);
                //std::cout<<"HERE"<<std::endl;
                Nor.resize(Ib.rows(), Fb.cols());
                compute_normals(Nor, P.transpose()*q2+x0, Ib, Fb);
                //std::cout<<Nor.rows()<<" "<<Ib.rows()<<std::endl;
                //std::cout<<Nor.cols()<<" "<<Ib.cols()<<std::endl;
                for(int i = 0; i < Ib.rows(); i++){
                    //f.segment<3>(3 * Ib(i)) -= Eigen::Vector3d(0.0, 10000.0, 0.0);
                    //std::cout<<"Inside loop"<<std::endl;
                    Eigen::RowVector3d n = Nor.row(i);
                    //std::cout<<"After n"<<std::endl;
                    //f.segment<3>(3 * Ib(i)) += c * n.transpose();
                    if(!n.hasNaN()){
                        //std::cout<<i<<std::endl;
                        Eigen::Vector3d H_tot;
                        Eigen::Vector3d p = (P.transpose()*q2+x0).segment<3>(3 * Ib(i));
                        if(bunny){
                            if(constant){
                                //std::cout<<"Inside constant"<<std::endl;
                                H << 0.0, 5000000.0, 0.0;// bunny
                            }
                            else{
                                bar_magnet(H, Po, p, 1.0 * 5e11);
                            }
                        }
                        else if(cube86){
                            if(constant){
                                H << 0.0, 250000.0, 0.0;
                                //H << 250000.0, 0.0, 0.0;
                            }
                            else{
                                bar_magnet(H, Po, p, 0.5 * 1e6);
                            }
                        }
                        else{
                            if(constant){
                                H << 0.0, 10000.0, 0.0;// arma
                            }
                            else{
                                bar_magnet(H, Po, p, 1.0 * 3e4);
                            }
                        }
                        H_tot = H - boundary_value.segment<3>(3 * i);
                        //H_tot = -boundary_value.segment<3>(3 * i);
                        c = (0.5 * mew * k * H_tot.transpose() * H_tot);
                        f.segment<3>(3 * Ib(i)) += c * n.transpose();
                    }
                }
                //std::cout<<"After loop"<<std::endl;
            }
            if(visualization){
                Visualize::add_scalar_field_visualization(f, fixed_point_indices);
            }
            f = P*f; 
        };

        //assemble stiffness matrix,
        auto stiffness = [&](Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q2, Eigen::Ref<const Eigen::VectorXd> qdot2) { 
            assemble_stiffness(K, P.transpose()*q2+x0, P.transpose()*qdot2, V, T, v0, C, D);
            K = P*K*P.transpose();
        };

        if(fully_implicit)
            implicit_euler(q, qdot, dt, M, energy, force, stiffness, tmp_qdot, tmp_force, tmp_stiffness);
        else
            linearly_implicit_euler(q, qdot, dt, M, force, stiffness, tmp_force, tmp_stiffness);
        

        
    KE = 0;
    PE = 0;

    for(unsigned int ei=0; ei<T.rows(); ++ei) {
        T_linear_tetrahedron(T_ele, P.transpose()*qdot, T.row(ei), density, v0(ei));
        KE += T_ele;

        V_linear_tetrahedron(V_ele, P.transpose()*q+x0, V, T.row(ei), v0(ei), C, D);
        PE += V_ele;
    }
    
    Visualize::add_energy(t, KE, PE);
}

inline void draw(Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, double t) {

    //update vertex positions using simulation
    Visualize::update_vertex_positions(0, P.transpose()*q + x0);

}

bool key_down_callback(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifiers) {

    if(key =='N') {
        std::cout<<"toggle integrators \n";
        fully_implicit = !fully_implicit;
    } else if(key == 'K') {
        
        skinning_on = !skinning_on;
        Visualize::toggle_skinning(skinning_on);
    }
    else if(key=='M'){
        magnet = ! magnet;
        std::cout<<"magnet = "<<magnet<<"\n";
    }

    else if(key=='V'){
        visualization = !visualization;
    }

    else if(key=='C'){
        constant = !constant;
    }

    else if(key=='Y'){
        Po(0, 1)+=mov;
        Visualize::update_point(Po);
    }

    else if(key=='H'){
        Po(0, 1)-=mov;
        Visualize::update_point(Po);
    }

    else if(key=='G'){
        Po(0, 0)-=mov;
        Visualize::update_point(Po);
    }

    else if(key=='J'){
        Po(0, 0)+=mov;
        Visualize::update_point(Po);
    }

    else if(key=='X'){
        Po(0, 2)-=mov;
        Visualize::update_point(Po);
    }

    else if(key=='B'){
        Po(0, 2)+=mov;
        Visualize::update_point(Po);
    }
    else if(key=='S'){
        menu = !menu;
    }
    //std::cout<<key<<std::endl;
    return false;
}
inline void assignment_setup(int argc, char **argv, Eigen::VectorXd &q, Eigen::VectorXd &qdot) {

    //load geometric data 

    if(argc > 1) {
        if(strcmp(argv[1], "arma") == 0) {
            read_tetgen(V,T, "../data/arma_6.node", "../data/arma_6.ele");
            igl::readOBJ("../data/armadillo.obj", V_skin, F_skin);
        
            bunny = false;
            fully_implicit = true;
            mov = 0.1;
        }
        else if(strcmp(argv[1], "cube86") == 0) {
            igl::readMESH("../data/cube86.mesh", V, T, F);

            cube86 = true;
            fully_implicit = true;
            mov = 0.1;
        }
        else if(strcmp(argv[1], "bunny") == 0){
            igl::readMESH("../data/coarser_bunny.mesh",V,T, F);
            igl::readOBJ("../data/bunny_skin.obj", V_skin, F_skin);

            bunny = true;
            fully_implicit = false;
            mov = 3.0;
        }
    }
    
    //std::cout<<V.rows()<<" "<<T.rows()<<std::endl;
    igl::boundary_facets(T, F);
    F = F.rowwise().reverse().eval();
    //std::cout<<F<<std::endl;
    Eigen::MatrixXi Ft = F.transpose();
    Eigen::VectorXi ind;
    ind = Eigen::Map<Eigen::VectorXi>(Ft.data(), Ft.rows()*Ft.cols());
    //std::cout<<ind<<std::endl;
    //std::set<int> o{ind.begin(), ind.end()};
    std::vector<int> tmp;
    for(int i = 0; i < ind.size(); i++) tmp.push_back(ind(i));
    std::sort(tmp.begin(), tmp.end());
    auto it = std::unique(tmp.begin(), tmp.end());
    tmp.resize(std::distance(tmp.begin(), it));

    //Eigen::VectorXi Ib((int)tmp.size());
    Ib.resize((int)tmp.size());
    for(int i = 0; i < tmp.size();i++) Ib(i) = tmp[i];
    //std::cout<<Ib.rows()<<std::endl;
    //std::cout<<Ib<<std::endl;

    Fb.resize(F.rows(), F.cols());

    for(int i = 0; i < F.rows(); i++)
    {
        for(int j = 0; j < F.cols(); j++)
        {
            int idx = F(i, j);
            int pos = -1;
            for(int k = 0; k < Ib.rows(); k++)
            {
                if(Ib(k) == idx)
                {
                    Fb(i, j) = k;
                    break;
                }
            }
            
        }
    }
    mi = V.colwise().minCoeff();
    mx = V.colwise().maxCoeff();
    double length = std::max(abs(mi(0) - mx(0)), std::max(abs(mi(1) - mx(1)), abs(mi(2) - mx(2))));
    length = round(length) + grid_length;
    cell_width = length / grid_length;
    //std::cout<<cell_width<<std::endl;
    //std::cout<<Fb<<std::endl;
    //std::cout<<Ib.rows()<<std::endl;
    //std::cout<<V.rows()<<std::endl;
    //for(int i = 0; i < tmp.size(); i++) std::cout<<tmp[i]<<" ";
    //std::cout<<std::endl;
    //auto it = std::unique(ind.begin(), ind.end());

    build_skinning_matrix(N, V, T, V_skin);

    //setup simulation 
    init_state(q,qdot,V);
    Eigen::MatrixXd Nor;
    Nor.resize(Ib.rows(), Fb.cols());
    compute_normals(Nor, q, Ib, Fb);
    //std::cout<<Nor.rows()<<" "<<Ib.rows()<<std::endl;
    //std::cout<<q.segment<3>(498*3)<<std::endl;
    //return;
    /*Eigen::MatrixXd Nor;
    compute_normals(Nor, q, Ib, Fb);
    std::cout<<Nor<<std::endl;*/

    //add geometry to scene
    Visualize::add_object_to_scene(V, F, V_skin, F_skin, N, Eigen::RowVector3d(244,165,130)/255.);
    Visualize::toggle_skinning(false);
    if(bunny){
        Po<<mx(0), mx(1), mx(2);
        Ini<<mx(0), mx(1), mx(2);
    }
    else if(cube86){
        Po<<mi(0) - 1.0, mi(1) + 1.5, mi(2) + 0.5;
        Ini<<mi(0) - 1.0, mi(1) + 1.5, mi(2) + 0.5;
        //Po<<mi(0) , mi(1), mi(2);
    }
    else{
        Po<<mi(0) + 0.45, mx(1) + 0.3, mi(2) + 0.2;
        Ini<<mi(0) + 0.45, mx(1) + 0.3, mi(2) + 0.2;
    }
    Visualize::add_point(Po);

    //Visualize::add_object_to_scene(Po, Pf, Pv_skin, Pf_skin, Pn, Eigen::RowVector3d(255,0,0)/255.);
    /*Eigen::MatrixXd V_box(8,3);
    V_box <<
    mi(0), mi(1), mi(2),
    mx(0), mi(1), mi(2),
    mx(0), mx(1), mi(2),
    mi(0), mx(1), mi(2),
    mi(0), mi(1), mx(2),
    mx(0), mi(1), mx(2),
    mx(0), mx(1), mx(2),
    mi(0), mx(1), mx(2);

    // Edges of the bounding box
    Eigen::MatrixXi E_box(12,2);
    E_box <<
    0, 1,
    1, 2,
    2, 3,
    3, 0,
    4, 5,
    5, 6,
    6, 7,
    7, 4,
    0, 4,
    1, 5,
    2, 6,
    7 ,3;

    Visualize::add_boundary_box(V_box, E_box);*/
    
    //bunny
    if(bunny)
        Visualize::set_picking_tolerance(1.);
    else if(cube86)
        Visualize::set_picking_tolerance(1.);
    else
        Visualize::set_picking_tolerance(0.01);

    //volumes of all elements
    igl::volume(V,T, v0);

    //Mass Matrix
    mass_matrix_mesh(M, qdot, T, density, v0);
    
    if(M.rows() == 0) {
        std::cout<<"Mass Matrix not implemented, quitting \n";
        std::exit(0);
    }
    
    //setup constraint matrix
    if(bunny)
        find_min_vertices(fixed_point_indices, V, 3);
    else if(cube86){
        find_min_vertices(fixed_point_indices, V, 0.001);}
        //std::cout<<fixed_point_indices.size()<<std::endl;
    else
        find_min_vertices(fixed_point_indices, V, 0.1);
    //material properties
    //bunny
    if(bunny) {
        YM = 6e6; //young's modulus
        mu = 0.4; //poissons ratio
        D = 0.5*(YM*mu)/((1.0+mu)*(1.0-2.0*mu));
        C = 0.5*YM/(2.0*(1.0+mu));
        k_selected = 1e8;
    } 
    else if(cube86){
        YM = 6e5; //young's modulus
        mu = 0.4; //poissons ratio
        D = 0.5*(YM*mu)/((1.0+mu)*(1.0-2.0*mu));
        C = 0.5*YM/(2.0*(1.0+mu));
        k_selected = 1e5;
    }
    else {
        //arma
        YM = 6e5; //young's modulus
        mu = 0.4; //poissons ratio
        D = 0.5*(YM*mu)/((1.0+mu)*(1.0-2.0*mu));
        C = 0.5*YM/(2.0*(1.0+mu));
        k_selected = 1e5;
    }

    P.resize(q.rows(),q.rows());
    P.setIdentity();
    fixed_point_constraints(P, q.rows(), fixed_point_indices);
    
    x0 = q - P.transpose()*P*q; //vector x0 contains position of all fixed nodes, zero for everything else    
    //correct M, q and qdot so they are the right size
    q = P*q;
    qdot = P*qdot;
    M = P*M*P.transpose();
    //std::cout<<x0.rows()<<std::endl;
    //igl additional menu setup
    // Add content to the default menu window
    Visualize::viewer_menu().callback_draw_custom_window = [&]()
    {
        // Define next window position + size
        ImGui::SetNextWindowPos(ImVec2(180.f * Visualize::viewer_menu().menu_scaling(), 10), ImGuiSetCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(800, 500), ImGuiSetCond_FirstUseEver);
        ImGui::Begin(
            "Energy Plot", nullptr,
            ImGuiWindowFlags_NoSavedSettings

        );

        ImVec2 min = ImGui::GetWindowContentRegionMin();
        ImVec2 max = ImGui::GetWindowContentRegionMax();

        max.x = ( max.x - min.x ) / 2;
        max.y -= min.y + ImGui::GetItemsLineHeightWithSpacing() * 3;

        Visualize::plot_energy("T", 1, ImVec2(-15,10), ImVec2(0,2e8), ImGui::GetColorU32(ImGuiCol_PlotLines));
        Visualize::plot_energy("V", 2, ImVec2(-15,10), ImVec2(0,2e7), ImGui::GetColorU32(ImGuiCol_HeaderActive));
        Visualize::plot_energy("T+V", 3, ImVec2(-15,10), ImVec2(0,4e8), ImGui::GetColorU32(ImGuiCol_ColumnActive));

        ImGui::End();
    };

    Visualize::viewer().callback_key_down = key_down_callback;

}

#endif

