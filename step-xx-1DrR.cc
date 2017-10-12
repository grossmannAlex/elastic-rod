/* ---------------------------------------------------------------------
 *
 * Copyright (C) 
 *
 * Author: Alexander Groﬂmann, University of Erlangen-Nuernberg, 2017
 */



#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
// #include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/grid/manifold_lib.h>


#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>

#include <fstream>
#include <iostream>
#include <stdlib.h>     /* abs */
#include <math.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/base/timer.h>

// my header
#include "helper_functions.h"


namespace StepXX
{
	
  using namespace dealii;

 
template <int dim, int spacedim>
class ElasticProblem
{
	public:
	ElasticProblem ();
	~ElasticProblem ();
	void run ();

	private:
	void third_grid(
		const double 	&len,
		const double 	&radius,
		const int 		&n_of_ref);
	void setup_system ( const bool first_step );
	void assemble_system ();
	void solve ();
	void refine_grid ();
	void output_results (const unsigned int cycle) const;


	InnerEnergyFunction material_law;
	
	Triangulation<dim, spacedim>   triangulation;
	DoFHandler<dim, spacedim>      dof_handler;
	FESystem<dim, spacedim>        fe;

	ConstraintMatrix 		hanging_node_constraints;
	SparsityPattern      	sparsity_pattern;
	SparseMatrix<double> 	system_matrix;

	Vector<double> solution;
	Vector<double> solution_old;
	Vector<double> system_rhs;
};


template <int dim, int spacedim>
ElasticProblem<dim, spacedim>::ElasticProblem ()
	:
	dof_handler (triangulation),
	fe ( 	FE_Q<dim, spacedim>(1), 3,
			FE_Q<dim, spacedim>(1), 3
	)
	{}




template <int dim, int spacedim>
ElasticProblem<dim, spacedim>::~ElasticProblem ()
{
    dof_handler.clear ();
}


template <int dim, int spacedim>
void ElasticProblem<dim, spacedim>::third_grid(
	  const double &len,
	  const double &radius,
	  const int &n_of_ref )
{

	std::cout 	<< "    length of the rod: " << len 	<< std::endl
				<< "    radius of the rod: " << radius 	<< std::endl
				<< std::endl<< std::endl;
  
	GridGenerator::hyper_cube (triangulation);
// 	GridGenerator::cylinder ( triangulation, radius, len );
// 	GridTools::rotate ( -M_PI / 2 , 1, triangulation );
// 	const CylindricalManifold<3,3> cyl_manifold (
// 		2,
// 		1e-10
// 		);
// 
//     triangulation.set_manifold (0, cyl_manifold);
//     triangulation.set_all_manifold_ids(0);
    triangulation.refine_global (n_of_ref);
//     triangulation.set_manifold (0);
    
}
  

  
template <int dim, int spacedim>
void ElasticProblem<dim, spacedim>::setup_system (const bool first_step)
{
	if ( first_step ) 
	{
		dof_handler.distribute_dofs (fe);
		solution_old.reinit ( dof_handler.n_dofs() );
	}
	
	hanging_node_constraints.clear ();
	DoFTools::make_hanging_node_constraints (
		dof_handler,
        hanging_node_constraints);
	hanging_node_constraints.close ();
	
	solution.reinit (dof_handler.n_dofs());
    system_rhs.reinit (dof_handler.n_dofs());
	
	DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp );
	hanging_node_constraints.condense (dsp);
	sparsity_pattern.copy_from (dsp);
    system_matrix.reinit (sparsity_pattern);
}



template <int dim, int spacedim>
void ElasticProblem<dim, spacedim>::assemble_system ()
{
	/*
	 * Describing which boundarys are prescirbed
	 * r1,r2,r3, theta1,theta2, theta3 - most left
	 * r1,r2,r3, theta1,theta2, theta3 - most right
	 */
	const bool prescribed_boundary_neumann[2][6] =
	{
		false, false, false, false, false, false,
		false, false, true, true, true, true
// 		true, true, true, true, true, true 	
	};
	
// 	const FEValuesExtractors::Vector centerline (0);
// 	const FEValuesExtractors::Vector rotations (spacedim);
	const FEValuesExtractors::Scalar r1 (0);
	const FEValuesExtractors::Scalar r2 (1);
	const FEValuesExtractors::Scalar r3 (2);
	const FEValuesExtractors::Scalar theta1 (3);
	const FEValuesExtractors::Scalar theta2 (4);
	const FEValuesExtractors::Scalar theta3 (5);
	

	QGauss<dim>  quadrature_formula(2);
	FEValues<dim, spacedim> fe_values (fe, quadrature_formula,
		update_values |
		update_gradients |
		update_quadrature_points |
		update_JxW_values );
	
	QGauss<dim-1>  face_quadrature_formula(1);
	FEFaceValues<dim> fe_face_values (fe, face_quadrature_formula,
		update_values | 
		update_gradients |
		update_quadrature_points |
// 		update_normal_vectors |
		update_JxW_values
	);

	
	const unsigned int 	dofs_per_cell 	= fe.dofs_per_cell;
	const unsigned int 	n_q_points    	= quadrature_formula.size();
	const unsigned int 	n_face_q_points = face_quadrature_formula.size();
	

	FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
	Vector<double>       cell_rhs (dofs_per_cell);

	std::vector<types::global_dof_index> 	local_dof_indices ( dofs_per_cell );
	std::vector<Tensor<1, spacedim> > 		old_solution_grad ( n_q_points );
	
	
	// helpers
	Tensor<1,3> temp_a;
	Tensor<1,3> temp_b;
	Tensor<1,3> temp_c;
	
	
	typename DoFHandler<dim, spacedim>::active_cell_iterator 
		cell = dof_handler.begin_active(),
		endc = dof_handler.end();
    for (; cell!=endc; ++cell) // ELEMENT
      {
        cell_matrix = 0;
        cell_rhs 	= 0;
        fe_values.reinit (cell);
		
		// GAx, GAy, EA, EIx, EIy, GJ
		material_law.set_coefficients(
			1.0, 1.0, 1.0, 1.0, 1.0, 1.0
		);
		
		
		for (unsigned int q_point=0; q_point<n_q_points; ++q_point) // QUADRATUR 
		{
			// old solution
			const double s_qpoint = fe_values.quadrature_point (q_point)[0];
			Tensor <1,3> r_old ( {0,0, s_qpoint} );
			Tensor <1,3> dr_ds_old ( { 0, 0, 1 } );
			Tensor <1,3> theta_old ( { 0, 0, 0 } );
			Tensor <1,3> dtheta_ds_old ( { 0, 0, 0 } );
			
			
			material_law.set_solution(
				r_old, dr_ds_old, theta_old , dtheta_ds_old );
			
			
			// forces and moments
 			const Tensor <1,3> n_vec = material_law.get_n();
 			const Tensor <1,3> m_vec = material_law.get_m();
			
			
			for (unsigned int i=0; i<dofs_per_cell; ++i) // DOFS(i)
			{
				const Tensor<1,3> phi_ir_value( {
					fe_values[r1].value(i,q_point),
					fe_values[r2].value(i,q_point),
					fe_values[r3].value(i,q_point) 
				});
				
				const Tensor<1,3> phi_ir_grad( {
					fe_values[r1].gradient(i,q_point)[0],
					fe_values[r2].gradient(i,q_point)[0],
					fe_values[r3].gradient(i,q_point)[0] 
				});
				
				const Tensor<1,3> psi_iR_value( {
					fe_values[theta1].value(i,q_point),
					fe_values[theta2].value(i,q_point),
					fe_values[theta3].value(i,q_point) 
				});
				
				const Tensor<1,3> psi_iR_grad( {
					fe_values[theta1].gradient(i,q_point)[0],
					fe_values[theta2].gradient(i,q_point)[0],
					fe_values[theta3].gradient(i,q_point)[0] 
				});


		        // RHS: Force and moment densities
				cell_rhs(i) += ( (
						  phi_ir_value * f_n_hat ( s_qpoint ) 
						+ psi_iR_value * f_m_hat ( s_qpoint )
					) * fe_values.JxW( q_point ) );
				// --------------


				for (unsigned int j=0; j<dofs_per_cell; ++j) // DOFS(j)
				{
					const Tensor<1,3> phi_jr_grad( {
						fe_values[r1].gradient(j,q_point)[0],
						fe_values[r2].gradient(j,q_point)[0],
						fe_values[r3].gradient(j,q_point)[0] 
					});
					
					const Tensor<1,3> psi_jR_value( {
						fe_values[theta1].value(j,q_point),
						fe_values[theta2].value(j,q_point),
						fe_values[theta3].value(j,q_point) 
					});
					
					const Tensor<1,3> psi_jR_grad( {
						fe_values[theta1].gradient(j,q_point)[0],
						fe_values[theta2].gradient(j,q_point)[0],
						fe_values[theta3].gradient(j,q_point)[0] 
					});
					

					temp_a = cross_product_3d ( dr_ds_old, psi_jR_value );
					temp_a += phi_jr_grad;
					
					// 1st part:
					temp_b = cross_product_3d ( n_vec , psi_jR_value );
					temp_b *= -1;
					temp_c = material_law.get_R_Cmn_RT(1,1) * temp_a;
					temp_b += temp_c;
					temp_c = material_law.get_R_Cmn_RT(1,2) * psi_jR_grad;
					temp_b += temp_c;
					
					cell_matrix(i,j) += (
						( phi_ir_grad * temp_b ) * fe_values.JxW(q_point) );

					
					// 2nd part
					temp_c = cross_product_3d ( dr_ds_old , temp_b );
					temp_b = cross_product_3d ( n_vec , phi_jr_grad );
					
					cell_matrix(i,j) += ( 
						( psi_iR_value * ( temp_b - temp_c ) ) * fe_values.JxW(q_point) );
					
					
					// 3rd part 
					temp_b = cross_product_3d ( m_vec , psi_jR_value );
					temp_b *= -1;
					temp_c = material_law.get_R_Cmn_RT(2,1) * temp_a;
					temp_b += temp_c;
					temp_c = material_law.get_R_Cmn_RT(2,2) * psi_jR_grad;
					temp_b += temp_c;
					
					
					cell_matrix(i,j) += (
						( psi_iR_grad * temp_b ) * fe_values.JxW(q_point) );
					
					
				} // end j
			} // end i 
		} // end q

		
		// RHS: Neumann
		for (unsigned int face_number=0; face_number < GeometryInfo<dim>::faces_per_cell; ++face_number)
			if (cell->face(face_number)->at_boundary() ) 	// && (cell->face(face_number)->boundary_id() == my_neumann_id)
			{	
				/*
				 * cell->face(i)->at_boundary() is always 1 for faces at the boundary
				 * face_number is here 0 (most left) or 1 (most right)
				 */ 
				
				fe_face_values.reinit (cell, face_number);
				for (unsigned int q_point=0; q_point<n_face_q_points; ++q_point) // here mostly one
				{
					const double face_s_qpoint = fe_face_values.quadrature_point ( q_point )[0]; // s
					
					Tensor <1,3> r_old ( {0,0, face_s_qpoint} );
					Tensor <1,3> dr_ds_old ( { 0, 0, 1 } );
					Tensor <1,3> theta_old ( { 0, 0, 0 } );
					Tensor <1,3> dtheta_ds_old ( { 0, 0, 0 } );
					
					
					material_law.set_solution(
						r_old, dr_ds_old, theta_old , dtheta_ds_old );
			
					
					// forces and moments
					const Tensor <1,3> n_vec = material_law.get_n();
					const Tensor <1,3> m_vec = material_law.get_m();
			
					
					for (unsigned int i=0; i<dofs_per_cell; ++i)
					{
						const Tensor <1,3> face_phi_ir_value( {
							fe_face_values[r1].value(i, q_point),
							fe_face_values[r2].value(i, q_point), 
							fe_face_values[r3].value(i, q_point), 
						});
						
						const Tensor <1,3> face_psi_iR_value( {
							fe_face_values[theta1].value(i, q_point),
							fe_face_values[theta2].value(i, q_point), 
							fe_face_values[theta3].value(i, q_point), 
						});
						
						
						if ( prescribed_boundary_neumann[face_number][i] == true )
						{	
							/*
							 * adding prescribed force terms 
							 */
							const Tensor <1,3> n_pre ( {0,0,0} );
							const Tensor <1,3> m_pre ( {0,0,0} );
							cell_rhs(i) += (
								  face_phi_ir_value * n_pre 
								+ face_psi_iR_value * m_pre
							) * fe_face_values.JxW( q_point );
						}
						else
						{	
							/* 
							 * adding unknowns to linear system 
							 */
							for (unsigned int j=0; j<dofs_per_cell; ++j)
							{
								
								const Tensor<1,3> face_phi_jr_grad( {
									fe_face_values[r1].gradient(j,q_point)[0],
									fe_face_values[r2].gradient(j,q_point)[0],
									fe_face_values[r3].gradient(j,q_point)[0] 
								});
								
								const Tensor<1,3> face_psi_jR_value( {
									fe_face_values[theta1].value(j,q_point),
									fe_face_values[theta2].value(j,q_point),
									fe_face_values[theta3].value(j,q_point) 
								});
								
								const Tensor<1,3> face_psi_jR_grad( {
									fe_face_values[theta1].gradient(j,q_point)[0],
									fe_face_values[theta2].gradient(j,q_point)[0],
									fe_face_values[theta3].gradient(j,q_point)[0] 
								});
								

								
								temp_a = cross_product_3d ( dr_ds_old, face_psi_jR_value );
								temp_a += face_phi_jr_grad;
								
								// n
								temp_b = cross_product_3d ( n_vec , face_psi_jR_value );
								temp_b *= -1;
								temp_c = material_law.get_R_Cmn_RT(1,1) * temp_a;
								temp_b += temp_c;
								temp_c = material_law.get_R_Cmn_RT(1,2)* face_psi_jR_grad;
								temp_b += temp_c;
								
								cell_matrix(i,j) += (
									( face_phi_ir_value * temp_b ) * fe_face_values.JxW( q_point ) );

								// m  
								temp_b = cross_product_3d ( m_vec , face_psi_jR_value );
								temp_b *= -1;
								temp_c = material_law.get_R_Cmn_RT(2,1)* temp_a;
								temp_b += temp_c;
								temp_c = material_law.get_R_Cmn_RT(2,2) * face_psi_jR_grad;
								temp_b += temp_c;

								cell_matrix(i,j) += (
									( face_psi_iR_value * temp_b ) * fe_face_values.JxW( q_point ) );
								
							} // j
						} // if: rhs known - rhs unknown
					} // i  
				} // q_point  
			} // for faces - if boundary 


          
        cell->get_dof_indices (local_dof_indices);
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              system_matrix.add (local_dof_indices[i],
                                 local_dof_indices[j],
                                 cell_matrix(i,j));

            system_rhs(local_dof_indices[i]) += cell_rhs(i);
          }
	} // end element 


	hanging_node_constraints.condense (system_matrix);
	hanging_node_constraints.condense (system_rhs);

	// Idee: alle Masken vorbereiten und Boundary IDs indivduell setzen
	// 12 Optionen (links, rechts, 6 Komponenten)
	std::map<types::global_dof_index,double> boundary_values;
	FEValuesExtractors::Scalar used_mask[]= {r1, r2, r3, theta1, theta2, theta3 };
	for (int i=0; i< 6 ; i++)
	{
		// most left boundary
		if (prescribed_boundary_neumann[0][i] == false)
			VectorTools::interpolate_boundary_values (
				dof_handler, 0,
				ZeroFunction<spacedim>( 6 ),
				boundary_values,
				fe.component_mask( used_mask[i] )
			);
		
		// most right boundary
		if (prescribed_boundary_neumann[1][i] == false)
			VectorTools::interpolate_boundary_values (
				dof_handler, 1,
				ZeroFunction<spacedim>( 6 ),
				boundary_values,
				fe.component_mask( used_mask[i] )
			);
	}
	
	MatrixTools::apply_boundary_values (
		boundary_values, system_matrix, solution, system_rhs);
}




template <int dim, int spacedim>
void ElasticProblem<dim, spacedim>::solve ()
{
	SolverControl           solver_control (10000, 1e-12);
	SolverCG<>              cg (solver_control);

//     PreconditionSSOR<> preconditioner;
//     preconditioner.initialize(system_matrix, 1.2);
// 
//      cg.solve (system_matrix, solution, system_rhs,
//                preconditioner);
	

	std::ofstream of("sys_mat.txt");
	system_matrix.print_formatted ( of );
	
	
// 	// CG
// 	cg.solve (	system_matrix,
// 				solution,
// 				system_rhs,
// 				PreconditionIdentity() );
	
	
	std::cout  << "    Solving linear system... ";
	Timer timer;
	timer.start ();
	SparseDirectUMFPACK  A_direct;
	A_direct.initialize(system_matrix);
	A_direct.vmult (solution, system_rhs);
	timer.stop ();
	std::cout 	<< "done (" << timer () << "s)"
				<< std::endl << std::endl;
	
				
	hanging_node_constraints.distribute (solution);
}

  
/**
 * Output:
 * 1) only components
 * 2) vector values (does not work)
 */
  template <int dim, int spacedim>
  void ElasticProblem<dim, spacedim>::output_results (const unsigned int cycle) const
  {

	std::cout << " size of U: " << solution.size() << std::endl;
	  
	std::string filename = "solution-";
    filename += ('0' + cycle);
    Assert (cycle < 10, ExcInternalError());

    filename += ".vtk";
    std::ofstream output (filename.c_str());

	std::vector<std::string> solution_names;
	solution_names.push_back ("rx");
    solution_names.push_back ("ry");
 	solution_names.push_back ("rz"); 
	solution_names.push_back ("theta1");
    solution_names.push_back ("theta2");
 	solution_names.push_back ("theta3"); 
	
    DataOut<dim, DoFHandler<dim,spacedim>> data_out;
    data_out.attach_dof_handler (dof_handler);
	data_out.add_data_vector (solution, solution_names);
    data_out.build_patches ();
    data_out.write_vtk (output);
	
	
/*
 *  	This does not work!
 *  	Could not find out how to tell the DataOut object, that the first three scalars are a vector.
 * 		Guessing this has s.t. to do with the spacedim.
 * 
 * 2017-09-26: After a short discussion with some colleges i decided to stay with the component 
 * expression and use the calculation tool of paraview to get a useful output.
 * 
 * */
	{
	/*
 
	std::ostringstream filename_vec;
	filename_vec 	<< "vector-solution"
					<< ".vtk";
    std::ofstream output_vec (filename_vec.str().c_str());	
	
	std::vector<std::string> solution_names_vec;
	solution_names_vec.push_back("centerline");
// 	solution_names_vec.push_back("centerline");
// 	solution_names_vec.push_back("centerline");
	solution_names_vec.push_back("rotation_1");
	solution_names_vec.push_back("rotation_2");
	solution_names_vec.push_back("rotation_3");
	
	std::vector<DataComponentInterpretation::DataComponentInterpretation> data_component_interpretation
	( 3, DataComponentInterpretation::component_is_part_of_vector );
	data_component_interpretation.push_back ( DataComponentInterpretation::component_is_scalar );
	data_component_interpretation.push_back ( DataComponentInterpretation::component_is_scalar );
	data_component_interpretation.push_back ( DataComponentInterpretation::component_is_scalar );

	DataOut<dim, DoFHandler<dim>> data_out_vec;
    data_out_vec.attach_dof_handler (dof_handler);
    data_out_vec.add_data_vector (
			solution, solution_names_vec,
			DataOut<dim, DoFHandler<dim>>::type_dof_data,
			data_component_interpretation);
	data_out_vec.build_patches ();
    data_out_vec.write_vtk (output_vec);*/
  }
  }




  template <int dim, int spacedim>
  void ElasticProblem<dim, spacedim>::run ()
  {
        third_grid ( 1 , 0.2 , 5 );
	  
        std::cout << "   Number of active cells:       "
                  << triangulation.n_active_cells()
                  << std::endl << std::endl;

        setup_system ( true );

        std::cout 	<< "   Total number of degrees of freedom: "
					<< dof_handler.n_dofs() << std::endl				  
					<< "   Number of degrees of freedom per cell: "
					<< fe.dofs_per_cell << std::endl<< std::endl;

        assemble_system ();
        solve ();
        output_results (0);
  }
  
} // end namespace step8


int main ()
{
  try
    {
	
	test_InnerEnergy_function();
		
	StepXX::ElasticProblem<1,1> elastic_problem_3d;
	elastic_problem_3d.run ();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
