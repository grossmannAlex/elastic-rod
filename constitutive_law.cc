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
        })
{
    v_zero = transpose(R_zero) * r_prime_zero;
    k_zero = axial_vec(transpose(R_zero) * R_prime_zero);
}

constitutive_law::~constitutive_law() {
    // some content  
}

void constitutive_law::set_coefficients(double GAx, double GAy, double EA, double EIx, double EIy, double GJ) {
    double GIx = 0;
    double GIy = 0;
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
        const Tensor < 1, 3 > &r_centerline,
        const Tensor < 1, 3 > &drds_centerline,
        const Tensor < 1, 3 > &theta_angles,
        const Tensor < 1, 3 > &dthetads_angles
        ) {
    r = r_centerline;
    r_prime = drds_centerline;
    theta = theta_angles;
    theta_prime = dthetads_angles;

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

    return ( R * result);
}

Tensor < 1, 3 > constitutive_law::get_m() {
    Tensor < 1, 3 > result;
    result = get_Cmn(2, 1) * (v - v_zero);
    result += get_Cmn(2, 2) * (k - k_zero);

    return ( R * result);
}

void constitutive_law::print_member() {
    std::cout << "r: "
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
            << k << std::endl << std::endl;
}
