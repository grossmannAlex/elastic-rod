/* 
 * File:   constitutive_law.cc
 * Author: alex
 * 
 * Created on October 13, 2017, 2:25 PM
 */


#include "constitutive_law.h"

constitutive_law::constitutive_law() :
r_prime_zero({0, 0, 1}),
R_zero({
    {1, 0, 0},
    {0, 1, 0},
    {0, 0, 1}
}),
R_prime_zero({
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0}
}) {
    v_zero = transpose(R_zero) * r_prime_zero;
    k_zero = axial_vec(transpose(R_zero) * R_prime_zero);
}

constitutive_law::~constitutive_law() {
    // some content  
}

void constitutive_law::set_coefficients() {
    //    double GIx = 0;
    //    double GIy = 0;
    
    
    double E = 210000;
    double G = 70000;
    double A = 0.01/4 * M_PI;
    double I = pow(0.1,4) * M_PI/4;
    double J = 2 * I;
    double GAx = G * A;
    double GAy = G * A;
    double EA = E * A;
    double EIx = E * I;
    double EIy = E * I;
    double GJ = G * J;
    double EIxy = 0;

    C11[0][0] = GAx;
    C11[0][1] = 0;
    C11[0][2] = 0;
    C11[1][0] = 0;
    C11[1][1] = GAy;
    C11[1][2] = 0;
    C11[2][0] = 0;
    C11[2][1] = 0;
    C11[2][2] = EA;

    C12[0][0] = 0;
    C12[0][1] = 0;
    C12[0][2] = 0;
    C12[1][0] = 0;
    C12[1][1] = 0;
    C12[1][2] = 0;
    C12[2][0] = 0;
    C12[2][1] = 0;
    C12[2][2] = 0;

    C21[0][0] = 0;
    C21[0][1] = 0;
    C21[0][2] = 0;
    C21[1][0] = 0;
    C21[1][1] = 0;
    C21[1][2] = 0;
    C21[2][0] = 0;
    C21[2][1] = 0;
    C21[2][2] = 0;

    C22[0][0] = EIx;
    C22[0][1] = EIxy;
    C22[0][2] = 0;
    C22[1][0] = EIxy;
    C22[1][1] = EIy;
    C22[1][2] = 0;
    C22[2][0] = 0;
    C22[2][1] = 0;
    C22[2][2] = GJ;
}

void constitutive_law::set_solution(
        FEValues < 1, 1 > &fe_values,
        const unsigned int n_q_points
        ) {
    // extractors
    const FEValuesExtractors::Scalar r1(0);
    const FEValuesExtractors::Scalar r2(1);
    const FEValuesExtractors::Scalar r3(2);
    const FEValuesExtractors::Scalar theta1(3);
    const FEValuesExtractors::Scalar theta2(4);
    const FEValuesExtractors::Scalar theta3(5);

    old_solution_r1.resize(n_q_points);
    old_solution_r2.resize(n_q_points);
    old_solution_r3.resize(n_q_points);
    old_solution_theta1.resize(n_q_points);
    old_solution_theta2.resize(n_q_points);
    old_solution_theta3.resize(n_q_points);

    old_solution_dr1ds.resize(n_q_points);
    old_solution_dr2ds.resize(n_q_points);
    old_solution_dr3ds.resize(n_q_points);
    old_solution_dtheta1ds.resize(n_q_points);
    old_solution_dtheta2ds.resize(n_q_points);
    old_solution_dtheta3ds.resize(n_q_points);

    // old solution as current state
    fe_values[r1].get_function_values(solution_old, old_solution_r1);
    fe_values[r2].get_function_values(solution_old, old_solution_r2);
    fe_values[r3].get_function_values(solution_old, old_solution_r3);
    fe_values[theta1].get_function_values(solution_old, old_solution_theta1);
    fe_values[theta2].get_function_values(solution_old, old_solution_theta2);
    fe_values[theta3].get_function_values(solution_old, old_solution_theta3);

    fe_values[r1].get_function_gradients(solution_old, old_solution_dr1ds);
    fe_values[r2].get_function_gradients(solution_old, old_solution_dr2ds);
    fe_values[r3].get_function_gradients(solution_old, old_solution_dr3ds);
    fe_values[theta1].get_function_gradients(solution_old, old_solution_dtheta1ds);
    fe_values[theta2].get_function_gradients(solution_old, old_solution_dtheta2ds);
    fe_values[theta3].get_function_gradients(solution_old, old_solution_dtheta3ds);

}

void constitutive_law::set_old_solution(
        const Vector<double> &input ){

    solution_old = input;
}

void constitutive_law::set_solution(
        FEFaceValues < 1, 1 > &fe_values,
        const unsigned int n_q_points
        ) {
    // extractors
    const FEValuesExtractors::Scalar r1(0);
    const FEValuesExtractors::Scalar r2(1);
    const FEValuesExtractors::Scalar r3(2);
    const FEValuesExtractors::Scalar theta1(3);
    const FEValuesExtractors::Scalar theta2(4);
    const FEValuesExtractors::Scalar theta3(5);

    old_solution_r1.resize(n_q_points);
    old_solution_r2.resize(n_q_points);
    old_solution_r3.resize(n_q_points);
    old_solution_theta1.resize(n_q_points);
    old_solution_theta2.resize(n_q_points);
    old_solution_theta3.resize(n_q_points);

    old_solution_dr1ds.resize(n_q_points);
    old_solution_dr2ds.resize(n_q_points);
    old_solution_dr3ds.resize(n_q_points);
    old_solution_dtheta1ds.resize(n_q_points);
    old_solution_dtheta2ds.resize(n_q_points);
    old_solution_dtheta3ds.resize(n_q_points);

    // old solution as current state
    fe_values[r1].get_function_values(solution_old, old_solution_r1);
    fe_values[r2].get_function_values(solution_old, old_solution_r2);
    fe_values[r3].get_function_values(solution_old, old_solution_r3);
    fe_values[theta1].get_function_values(solution_old, old_solution_theta1);
    fe_values[theta2].get_function_values(solution_old, old_solution_theta2);
    fe_values[theta3].get_function_values(solution_old, old_solution_theta3);

    fe_values[r1].get_function_gradients(solution_old, old_solution_dr1ds);
    fe_values[r2].get_function_gradients(solution_old, old_solution_dr2ds);
    fe_values[r3].get_function_gradients(solution_old, old_solution_dr3ds);
    fe_values[theta1].get_function_gradients(solution_old, old_solution_dtheta1ds);
    fe_values[theta2].get_function_gradients(solution_old, old_solution_dtheta2ds);
    fe_values[theta3].get_function_gradients(solution_old, old_solution_dtheta3ds);

}

void constitutive_law::get_solution(
        const unsigned int q_point,
        const double s_qpoint
        ) {

    // get values of solution_old
    r[0] = old_solution_r1[q_point];
    r[1] = old_solution_r2[q_point];
    r[2] = old_solution_r3[q_point] + s_qpoint;

    r_prime[0] = old_solution_dr1ds[q_point][0];
    r_prime[1] = old_solution_dr2ds[q_point][0];
    r_prime[2] = old_solution_dr3ds[q_point][0] + 1.0;

    
    theta[0] = old_solution_theta1[q_point];
    theta[1] = old_solution_theta2[q_point];
    theta[2] = old_solution_theta3[q_point];

    theta_prime[0] = old_solution_dtheta1ds[q_point][0];
    theta_prime[1] = old_solution_dtheta2ds[q_point][0];
    theta_prime[2] = old_solution_dtheta3ds[q_point][0];

    Tensor < 2, 3 > K, R_prime;
    rodriguez_formula(theta, R);

    R_prime = R * skw_mat(theta_prime);

    v = transpose(R) * r_prime;
    K = transpose(R) * R_prime;
    k = axial_vec(K);
}

Tensor < 2, 3 > constitutive_law::get_Cmn(const int sec_m, const int sec_n) {
    Tensor < 2, 3 > result;

    switch (sec_m) {
        case 1:
            switch (sec_n) {
                case 1:
                    result = C11;
                    break;
                case 2:
                    result = C12;
                    break;
                default:
                    std::cout << "--- not found ---" << std::endl;
            }
            break;
        case 2:
            switch (sec_n) {
                case 1:
                    result = C21;
                    break;
                case 2:
                    result = C22;
                    break;
                default:
                    std::cout << "--- not found ---" << std::endl;
            }
            break;
        default:
            std::cout << "--- not found ---" << std::endl;
    }
    return result;
}

Tensor < 2, 3 > constitutive_law::get_R_Cmn_RT(const int sec_m, const int sec_n) {
    Tensor < 2, 3 > result;
    result = R * get_Cmn(sec_m, sec_n);
    result = result * transpose(R);
    return result;
}

Tensor < 1, 3 > constitutive_law::get_n() {
    Tensor < 1, 3 > result;
    result = get_Cmn(1, 1) * (v - v_zero);
    result += get_Cmn(1, 2) * (k - k_zero);
    result += get_Cmn(2, 1) * (k - k_zero);

    return ( R * result);
}

Tensor < 1, 3 > constitutive_law::get_m() {
    Tensor < 1, 3 > result;
    result = get_Cmn(2, 2) * (k - k_zero );
    result += get_Cmn(2, 1) * (v - v_zero);
    result += get_Cmn(1, 2) * (v - v_zero);

    return ( R * result);
}

void constitutive_law::print_member() {
    std::cout << std::endl
            << "r: "
            << r << std::endl
            << "r': "
            << r_prime << std::endl
            << "theta: "
            << theta << std::endl
            << "theta': "
            << theta_prime << std::endl
            << "R: "
            << R << std::endl
            << "v: "
            << v << std::endl
            << "k: "
            << k << std::endl 
            << "n: "
            << get_n() << std::endl
            << "m: "
            << get_m() << std::endl
            << "psi: "
            << psi() << std::endl << std::endl;
            
}

double constitutive_law::psi() {
    double result = 0;
    
    Tensor <1,3> temp;
    Tensor <1,3> v_star = v - v_zero;
    
    temp = get_Cmn(1,1) * v_star;
    result += 0.5 * ( v_star * temp );
    
    temp = get_Cmn(2,2) * k;
    result += 0.5 * ( k * temp );
    
    temp = get_Cmn(1,2) * v_star;
    result += 0.5 * ( k * temp );
    
    temp = get_Cmn(2,1) * k;
    result += 0.5 * ( v_star * temp );
    
    return result;
}
