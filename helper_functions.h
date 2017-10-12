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

#define PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647093844609550582231725359408128481117450284102701938521105559644622948954930381964428810975665933446128475648233786783165271201909145648566923460348610454326648213393607260249141273

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
		 { 0, 1, 0 });
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
Tensor<1,3> axial_vec( const Tensor<2,3> &inp )
{
	Tensor <1,3> vec;
	vec[0] = inp[2][1];
	vec[1] = inp[0][2];
	vec[2] = inp[1][0];
	
	// test
	Tensor <2,3> test;
	test = inp + transpose(inp);
	if ( test.norm() > 0 )
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
	 
	if ( inp_vector.norm() == 0 ) 
		{ return; }
	 
	 
	Tensor<2,3> skw_matrix;
	skw_matrix = skw_mat( inp_vector );
	
	const long double angle = inp_vector.norm();
	
	out_matrix += ( 
		( sin( angle ) * skw_matrix )
		/ angle
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

	/* set solution of newton iteration oder inital configuration */
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
	
	
	/* print private members */
	void print_member();
	
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
		
		const Tensor <1,3> r_prime_zero;
		const Tensor <2,3> R_zero;
		const Tensor <2,3> R_prime_zero;
		Tensor <1,3> v_zero;
		Tensor <1,3> k_zero;
};

  


InnerEnergyFunction::InnerEnergyFunction() :
	r_prime_zero( {0,0,1} ), 
	R_zero( { {1,0,0}, {0,1,0}, {0,0,1} } ),
	R_prime_zero( { {0,0,0}, {0,0,0}, {0,0,0} } )
{
	v_zero = transpose(R_zero) * r_prime_zero; 
	k_zero = axial_vec ( transpose(R_zero) * R_prime_zero );
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
	result = get_Cmn(1,1) * (v - v_zero);
	result += get_Cmn(1,2) * (k - k_zero);
	
	return ( R * result );
}



Tensor< 1, 3 > InnerEnergyFunction::get_m()
{
	Tensor< 1, 3 > result;
	result = get_Cmn(2,1) * (v - v_zero);
	result += get_Cmn(2,2) * (k - k_zero);
	
	return ( R * result );
}



void InnerEnergyFunction::print_member()
{
		std::cout 	<< "r: " 		
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
					<< k << std::endl<< std::endl;
}




/*
 * some tests
 */
void test_InnerEnergy_function()
{
	std::cout << "running some unit-tests ... " << std::endl;
	
	InnerEnergyFunction my_func;
	Tensor <1,3> 
		r ({0,0,1}),
		r_prime ({0,0,1}), 
		theta ({0,0,PI/4}), 
		theta_prime ({0,0,0});
	
	my_func.set_solution( r, r_prime, theta, theta_prime );
	
	my_func.set_coefficients(
		1.0, 1.0, 1.0, 1.0, 1.0, 1.0
	);
	
	std::cout 	<< my_func.get_Cmn(1,1) << std::endl
				<< my_func.get_Cmn(1,2) << std::endl
				<< my_func.get_Cmn(2,1) << std::endl
				<< my_func.get_Cmn(2,2) << std::endl<<std::endl;
	
	std::cout 	<< my_func.get_R_Cmn_RT(1,1) << std::endl
				<< my_func.get_R_Cmn_RT(1,2) << std::endl
				<< my_func.get_R_Cmn_RT(2,1) << std::endl
				<< my_func.get_R_Cmn_RT(2,2) << std::endl<<std::endl;
				
	std::cout 	<< "n: " << my_func.get_n() << std::endl
				<< "m: " << my_func.get_m() << std::endl<< std::endl;
				
	my_func.print_member();
	
	
	std::cout << "... passing all of them." << std::endl << std::endl;
}
