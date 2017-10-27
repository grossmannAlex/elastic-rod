/* 
 * File:   distributed_loads.cc
 * Author: alex
 * 
 * Created on October 13, 2017, 2:20 PM
 */

#include "distributed_loads.h"

distributed_loads::distributed_loads()
{
}

//distributed_loads::distributed_loads(const distributed_loads& orig) {
//}

distributed_loads::~distributed_loads() {
}

Tensor < 1, 3 > distributed_loads::f_n_hat(
        const double &s,
        const Tensor < 2, 3 > R) {
    assert(s > 0);
    using namespace dealii;
    Tensor < 1, 3 > result;
    result[0] = 0.0;
    result[1] = 0.0;
    result[2] = 0.0;

    Tensor < 1, 3 > Rn_pre;
    Rn_pre[0] = 0.0;
    Rn_pre[1] = 0.0;
    Rn_pre[2] = 0.0;
    Rn_pre = R * Rn_pre;
    return (result + Rn_pre)*alpha;
}

Tensor < 1, 3 > distributed_loads::f_m_hat(
        const double &s,
        const Tensor < 2, 3 > R) {
    assert(s > 0);
    Tensor < 1, 3 > result;
    result[0] = 0.0;
    result[1] = 0.0;
    result[2] = 0.0;

    Tensor < 1, 3 > Rm_pre;
    Rm_pre[0] = 0.0;
    Rm_pre[1] = 0.0;
    Rm_pre[2] = 0.0;
    Rm_pre = R * Rm_pre;
    return (result + Rm_pre)*alpha;
}

Tensor < 1, 3 > distributed_loads::f_n_pre(
        const Tensor < 2, 3 > R) {
    Tensor < 1, 3 > n_pre;
    n_pre[0] = 0.0;
    n_pre[1] = 0.0;
    n_pre[2] = 0.0;

    Tensor < 1, 3 > Rn_pre;
    Rn_pre[0] = 400.0;
    Rn_pre[1] = 0.0;
    Rn_pre[2] = 0.0;
    Rn_pre = R * Rn_pre;
    return (n_pre + Rn_pre)*alpha;
}

Tensor < 1, 3 > distributed_loads::f_m_pre(
        const Tensor < 2, 3 > R) {
    Tensor < 1, 3 > m_pre;
    m_pre[0] = 0.0;
//    m_pre[1] = 1.25*M_PI*100;
    m_pre[1] = 0.0;
    m_pre[2] = 0.0;

    Tensor < 1, 3 > Rm_pre;
    Rm_pre[0] = 0.0;
    Rm_pre[1] = 0.0;
    Rm_pre[2] = 0.0;
    Rm_pre = R * Rm_pre;
    return (m_pre + Rm_pre)*alpha;
}


bool distributed_loads::ask_for_boundary(
        const unsigned int i,
        const unsigned int j) {
    
    bool prescribed_boundary_neumann[2][6];
    /*
     * Describing which boundarys are prescirbed
     * r1,r2,r3, theta1,theta2, theta3 - most left
     * r1,r2,r3, theta1,theta2, theta3 - most right
     */
    prescribed_boundary_neumann[0][0] = false;
    prescribed_boundary_neumann[0][1] = false;
    prescribed_boundary_neumann[0][2] = false;
    prescribed_boundary_neumann[0][3] = false;
    prescribed_boundary_neumann[0][4] = false;
    prescribed_boundary_neumann[0][5] = false;

//    prescribed_boundary_neumann[1][0] = false;
//    prescribed_boundary_neumann[1][1] = false;
//    prescribed_boundary_neumann[1][2] = false;
//    prescribed_boundary_neumann[1][3] = false;
//    prescribed_boundary_neumann[1][4] = false;
//    prescribed_boundary_neumann[1][5] = false;

    prescribed_boundary_neumann[1][0] = true;
    prescribed_boundary_neumann[1][1] = true;
    prescribed_boundary_neumann[1][2] = true;
    prescribed_boundary_neumann[1][3] = true;
    prescribed_boundary_neumann[1][4] = true;
    prescribed_boundary_neumann[1][5] = true;
    
    return prescribed_boundary_neumann[i][j];
}
    
double distributed_loads::get_dirichlet_left( int entry_num )
{
    double most_left_dirichlet[6];
    
    // most left boundary condition
    
//    std::string msg = "distributed_loads::get_dirichlet_left - error \n\tneumann and dirichlet boundary is requested for the same point";
//    std::ostream my_stream( msg.c_str() );
    
    std::string my_string = "neumann and dirichlet prescribed for the same boundary";
    Assert( ask_for_boundary(0, entry_num) == false , ExcMessage( my_string ) ); 
            
     
    most_left_dirichlet[0] = 0;
    most_left_dirichlet[1] = 0;
    most_left_dirichlet[2] = 0;
    most_left_dirichlet[3] = 0;
    most_left_dirichlet[4] = 0;
    most_left_dirichlet[5] = 0;
    
    return most_left_dirichlet[entry_num];
}

double distributed_loads::get_dirichlet_right( int entry_num )
{
    double most_right_dirichlet[6];
    
    // most left boundary condition
    std::string my_string = "neumann and dirichlet prescribed for the same boundary";
    Assert( ask_for_boundary(1, entry_num) == false , ExcMessage( my_string ) ); 
    
    most_right_dirichlet[0] = 0;
    most_right_dirichlet[1] = 0;
    most_right_dirichlet[2] = 0;
    most_right_dirichlet[3] = 0;
    most_right_dirichlet[4] = 0;
    most_right_dirichlet[5] = 0;
    
    return most_right_dirichlet[entry_num];
}


void distributed_loads::set_initial_num_steps( int num ){
    Assert( num > 0 , ExcMessage( " step number: must be positive " ) ); 
    initial_num_steps = num;
    delta_alpha = 1.0/num;
}

void distributed_loads::set_alpha_zero( double num ){
    Assert( num >= 0 , ExcMessage( " alpha zero: must be positive or zero " ) ); 
    alpha_zero = num;
}

void distributed_loads::decrease_step(){
    std::cout << "\n\n decrease load step \n\n";
    beta *= 2;
}

void distributed_loads::increase_step(){
    std::cout << "\n\n increase load step ";
    if (beta > 1)
        beta /= 2;
    else
        std::cout << "already max";
    
    std::cout << "\n\n";
}

void distributed_loads::increase_alpha(){
    alpha = alpha_zero + delta_alpha / beta;
    
    if (alpha > 1.0)
        alpha = 1.0;
    
    alpha_zero = alpha;
    
    std::cout
            << "\n\n increas load to "
            << alpha
            << std::endl << std::endl;
}

void distributed_loads::reset_beta(){
    std::cout << "\n reset beta \n";
    beta = 1;
}

double distributed_loads::get_alpha(){
    return alpha;
}
    
//
////bool distributed_loads::decrease_load()
////{
////    Assert( load_factor > 0 , ExcMessage( " no " ) ); 
////    load_factor *= 0.9;
////    
////    std::cout << "\ndecrease load\n";
////    return false;
////}
//
//bool distributed_loads::increase_load()
//{
////    Assert( load_factor > 0 , ExcMessage( " no " ) ); 
////   
////        load_factor /= 0.9;
////    
////    if ( load_factor > 1.0)
////    {
////        load_factor = 1.0;
////        std::cout << "\nmax load\n";
////        return false;
////    }
////    std::cout << "\nincrease load\n";
////    return true;
//    
//    if (load_steps > 1){
//        load_steps--;
//        std::cout << "\nincrease load\n";
//        load_factor = 1.0/load_steps;
//        return true;
//    }
//    std::cout << "\nincrease load - max\n";
//    return false;
//}
//
//void distributed_loads::set_load_steps(int num)
//{
//    load_steps = num+1;
//    load_factor = 0.0;
//}
//
//double distributed_loads::get_load_factor()
//{
//    return load_factor;
//}
//
