/*
 * <one line to give the program's name and a brief idea of what it does.>
 * Copyright (C) 2017  <copyright holder> <email>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

// #ifndef SKW_MATRIX_H
// #define SKW_MATRIX_H
// 
// class skw_matrix
// {
// public:
//     skw_matrix(  );
//     ~skw_matrix();
// };
// 
// 
// #endif // SKW_MATRIX_H


#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/diagonal_matrix.h>

#define PI 3.14159265

using namespace dealii;


// template <int dim, int spacedim>
// void right_hand_side (const std::vector<Point<dim> > &points,
//                         std::vector<Tensor<1, spacedim> >   &values)
// {
// 	for (unsigned int point_n = 0; point_n < points.size(); ++point_n)
// 	{
// 		values[point_n][0] = 5;
// 		values[point_n][1] = 1;
// 		values[point_n][2] = 1;
// 	}    
// }


// force density
Tensor<1,3> f_n_hat ( const double &s )
{
	 const Tensor<1,3> result(
		 { 100, 0, 0 });
	return result;
}


// moment density
Tensor<1,3> f_m_hat ( const double &s )
{
	 const Tensor<1,3> result(
		 { 0, 0, 0 });
	return result;
}


/*
 * building a skew symmetric matrix of a vector
 */
Tensor<2,3> skw_mat( const Tensor<1,3> &inp )
{
	Tensor <2,3> skw;
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
Tensor<1,3> axial_vec( const Tensor<2,3> &inp )
{
	Tensor <1,3> vec;
	vec[0] = inp[2][1];
	vec[1] = inp[0][2];
	vec[2] = inp[1][0];
	
	// test
	Tensor <2,3> test;
	test = inp + transpose(inp);
	if ( test.norm_square() > 0 )
		std::cout << "--- skew-symmetric matrix is not skew symmetric! ---" <<std::endl;
	
	return vec;
}



/*
 * Rodrigues Formula:
 * R = exp(A)
 * A = -A^T can be build of inp_vector=a
 * A = A(a)
 * exp(A) = I + sin( norm(a) ) A/norm(a) + (1-cos( norm(a) )) [ a/norm(a) ]^2
 */
void rodriguez_formula ( const Tensor<1,3> &inp_vector , Tensor<2,3> &out_matrix )
{
	out_matrix.clear();
	out_matrix[0][0] = 1;
	out_matrix[1][1] = 1;
	out_matrix[2][2] = 1;
	 
	if ( inp_vector.norm_square() == 0 ) 
		{ return; }
	 
	 
	Tensor<2,3> skw_matrix;
	skw_matrix = skw_mat( inp_vector );
	
	const double angle = inp_vector.norm_square();
	
	out_matrix += ( 
		( sin( angle ) * skw_matrix )
		/ inp_vector.norm_square()
	);
	
	out_matrix += ( 
		( ( 1.0 - cos( angle ) ) * skw_matrix * skw_matrix )
		/ ( angle * angle )
	);
}



/*
 * ....
 */
class InnerEnergyFunction
{
public:
	InnerEnergyFunction ();
	~InnerEnergyFunction ();

	
	/* set initial values ... */
	void set_solution ();
	/* set solution of newton iteration */
	void set_solution ( 
		const Tensor <1,3> &r_centerline ,
		const Tensor <1,3> &drds_centerline ,
		const Tensor <1,3> &theta_angles ,
		const Tensor <1,3> &dthetads_angles
	);
	
	
	/* material parameter */
	void set_coefficients 	(
		double GAx, double GAy, double EA, double EIx, double EIy, double GJ
	);
	
	
	/* choose a material law */
	void set_constitutive_law ( const std::string &name_of_law );
	
	
	/* R Matrix */
	Tensor <2,3> get_R() { return R; };
	
	
	/* constitutive law */
	Tensor <2,3> get_Cmn 		( const int sec_m, const int sec_n );
	Tensor <2,3> get_R_Cmn_RT 	( const int sec_m, const int sec_n );
	
	
	/* dependet quantities for current state */
	Tensor <1,3> get_n();
	Tensor <1,3> get_m();
	
	
	private:
		std::string constitutive_law;
		
		/* r,r', theta, theta' */
		Tensor <1,3> r;
		Tensor <1,3> r_prime;
		Tensor <1,3> theta;
		Tensor <1,3> theta_prime;
		Tensor <2,3> R;
		
		
		/* coefficients */
		Tensor <2,3> C11;
		Tensor <2,3> C12;
		Tensor <2,3> C21;
		Tensor <2,3> C22;
	
		Tensor <2,3> R_C11_RT;
		Tensor <2,3> R_C12_RT;
		Tensor <2,3> R_C21_RT;
		Tensor <2,3> R_C22_RT;
		
		Tensor <1,3> v; 	// shear strain
		Tensor <1,3> k; 	// curvature
};

  


InnerEnergyFunction::InnerEnergyFunction () 
{
	// some content 
}


InnerEnergyFunction::~InnerEnergyFunction ()
{
	// some content  
}


void InnerEnergyFunction::set_coefficients(double GAx, double GAy, double EA, double EIx, double EIy, double GJ)
{
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
	C12[0][2] = -GIx;
	C12[1][0] = 0;
	C12[1][1] = 0;
	C12[1][2] = GIy;
	C12[2][0] = EIx;
	C12[2][1] = -EIy;
	C12[2][2] = 0;
	
	C21[0][0] = 0;
	C21[0][1] = 0;
	C21[0][2] = EIx;
	C21[1][0] = 0;
	C21[1][1] = 0;
	C21[1][2] = -EIy;
	C21[2][0] = -GIx;
	C21[2][1] = GIy;
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


void InnerEnergyFunction::set_solution ( 
		const Tensor <1,3> &r_centerline ,
		const Tensor <1,3> &drds_centerline ,
		const Tensor <1,3> &theta_angles ,
		const Tensor <1,3> &dthetads_angles
	)
{
	r 			= r_centerline;
	r_prime 	= drds_centerline;
	theta		= theta_angles;
	theta_prime = dthetads_angles;
	
	Tensor<2,3> K, R_prime;
	rodriguez_formula( theta , R );
	
	R_prime = R * skw_mat( theta_prime );

	
	v = transpose(R) * r_prime;
	K = transpose(R) * R_prime;
	k = axial_vec( K );
}


Tensor< 2, 3 > InnerEnergyFunction::get_Cmn ( const int sec_m, const int sec_n )
{
	Tensor< 2, 3 > result;
	
	switch (sec_m)
	{
		case 1:
			switch (sec_n)
			{
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
			switch (sec_n)
			{
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



Tensor< 2, 3 > InnerEnergyFunction::get_R_Cmn_RT( const int sec_m, const int sec_n )
{
	Tensor< 2, 3 > result;
	result = R * get_Cmn( sec_m, sec_n );
	result = result * transpose(R);
	return result;
}



Tensor< 1, 3 > InnerEnergyFunction::get_n()
{
	Tensor< 1, 3 > result;
	result = get_Cmn(1,1) * v;
	result += get_Cmn(1,2) * k;
	
	return ( R * result );
}



Tensor< 1, 3 > InnerEnergyFunction::get_m()
{
	Tensor< 1, 3 > result;
	result = get_Cmn(2,1) * v;
	result += get_Cmn(2,2) * k;
	
	return ( R * result );
}

