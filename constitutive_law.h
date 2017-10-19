/* 
 * File:   constitutive_law.h
 * Author: alex
 *
 * Created on October 13, 2017, 2:25 PM
 */

#ifndef CONSTITUTIVE_LAW_H
#define	CONSTITUTIVE_LAW_H

#include <deal.II/lac/vector.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/diagonal_matrix.h>
#include <math.h>       /* pow */

using namespace dealii;

/*
 * building a skew symmetric matrix of a vector
 */
Tensor < 2, 3 > skw_mat(const Tensor < 1, 3 > &inp) {
    Tensor < 2, 3 > skw;
    skw[0][0] = 0.0;
    skw[1][1] = 0.0;
    skw[2][2] = 0.0;

    skw[0][1] = -inp[2];
    skw[0][2] = inp[1];
    skw[1][0] = inp[2];
    skw[1][2] = -inp[0];
    skw[2][0] = -inp[1];
    skw[2][1] = inp[0];

    return skw;
}

/*
 * building a axial vec of a skew-symmetric matrix
 */
Tensor < 1, 3 > axial_vec(const Tensor < 2, 3 > &inp) {
    Tensor < 1, 3 > vec;
    vec[0] = inp[2][1];
    vec[1] = inp[0][2];
    vec[2] = inp[1][0];

    // test
    Tensor < 2, 3 > test;
    test = inp + transpose(inp);
    if (test.norm() > 1.0e-12)
        std::cout
                << "--- skew-symmetric matrix is not skew symmetric! --- "
                << test.norm() << " ---"
                << std::endl;

    return vec;
}

/*
 * Rodrigues Formula:
 * R = exp(A)
 * A = -A^T can be build of inp_vector=a
 * A = A(a)
 * exp(A) = I + sin( norm(a) ) A/norm(a) + (1-cos( norm(a) )) [ a/norm(a) ]^2
 */
void rodriguez_formula(const Tensor < 1, 3 > &inp_vector, Tensor < 2, 3 > &out_matrix) {
    out_matrix.clear();
    out_matrix[0][0] = 1;
    out_matrix[1][1] = 1;
    out_matrix[2][2] = 1;

    if (inp_vector.norm() == 0) {
        return;
    }


    Tensor < 2, 3 > skw_matrix;
    skw_matrix = skw_mat(inp_vector);

    const long double angle = inp_vector.norm();

    out_matrix += (
            (sin(angle) * skw_matrix)
            / angle
            );

    out_matrix += (
            ((1.0 - cos(angle)) * skw_matrix * skw_matrix)
            / (angle * angle)
            );
}

class constitutive_law {
public:
    constitutive_law();
    ~constitutive_law();

    /* material parameter */
    void set_coefficients();

    void set_old_solution( const Vector<double> &input );

    /* set solution of newton iteration oder inital configuration */
    void set_solution(
            FEValues < 1, 1 > &fe_values,
            const unsigned int n_q_points
            );
    
    void set_solution(
            FEFaceValues < 1, 1 > &fe_values,
            const unsigned int n_q_points
            );

    void get_solution(
            const unsigned int q_point,
            const double s_qpoint
            );

    double psi();


    /* print private members */
    void print_member();

    /* R Matrix */
    Tensor < 2, 3 > get_R() {
        return R;
    };
    
    Tensor < 1, 3 > get_drds() {
        return r_prime;
    };


    /* constitutive law */
    Tensor < 2, 3 > get_Cmn(const int sec_m, const int sec_n);
    Tensor < 2, 3 > get_R_Cmn_RT(const int sec_m, const int sec_n);


    /* dependet quantities for current state */
    Tensor < 1, 3 > get_n();
    Tensor < 1, 3 > get_m();


private:

    /* r,r', theta, theta' */
    Tensor < 1, 3 > r;
    Tensor < 1, 3 > r_prime;
    Tensor < 1, 3 > theta;
    Tensor < 1, 3 > theta_prime;
    Tensor < 2, 3 > R;


    /* coefficients */
    Tensor < 2, 3 > C11;
    Tensor < 2, 3 > C12;
    Tensor < 2, 3 > C21;
    Tensor < 2, 3 > C22;

    Tensor < 2, 3 > R_C11_RT;
    Tensor < 2, 3 > R_C12_RT;
    Tensor < 2, 3 > R_C21_RT;
    Tensor < 2, 3 > R_C22_RT;

    Tensor < 1, 3 > v; // shear strain
    Tensor < 1, 3 > k; // curvature


    Vector<double> solution_old;
    
    std::vector<double> old_solution_r1;
    std::vector<double> old_solution_r2;
    std::vector<double> old_solution_r3;
    std::vector<double> old_solution_theta1;
    std::vector<double> old_solution_theta2;
    std::vector<double> old_solution_theta3;

    std::vector<Tensor < 1, 1 > > old_solution_dr1ds;
    std::vector<Tensor < 1, 1 > > old_solution_dr2ds;
    std::vector<Tensor < 1, 1 > > old_solution_dr3ds;
    std::vector<Tensor < 1, 1 > > old_solution_dtheta1ds;
    std::vector<Tensor < 1, 1 > > old_solution_dtheta2ds;
    std::vector<Tensor < 1, 1 > > old_solution_dtheta3ds;

    /* 
     * This is somehow special!
     * Imagine a rod like r(s,t=0) = [0,0,s]
     * therefore r'(s,t=0) is (0,0,1)
     * Assume in that cas the stored energy is zero,
     * so one would expect the strains also to be zero
     * by definition v = R^T r' = [0,0,1] (here)
     * but looking to the forces due to shear strains 
     * n = R * (dpsi / dv) should be [0, 0, 0]
     * the material law need therefore something like
     * r' - r'(0) to hold former expectation
     */

    const Tensor < 1, 3 > r_prime_zero;
    const Tensor < 2, 3 > R_zero;
    const Tensor < 2, 3 > R_prime_zero;
    Tensor < 1, 3 > v_zero;
    Tensor < 1, 3 > k_zero;
};




#endif	/* CONSTITUTIVE_LAW_H */

