/* 
 * File:   elastic_rod.cc
 * Author: alex
 * 
 * Created on October 13, 2017, 2:24 PM
 */

#include "elastic_rod.h"


using namespace dealii;

template <int dim, int spacedim>
elastic_rod<dim, spacedim>::elastic_rod()
:
dof_handler(triangulation),
fe(     FE_Q<dim, spacedim>(2), 3,
        FE_Q<dim, spacedim>(1), 3) ,
        quadrature(2),
        face_quadrature(1),
        max_newton_iter(500),
        global_refinement(4),
        global_alpha_off(false)
{
}

template <int dim, int spacedim>
elastic_rod<dim, spacedim>::~elastic_rod() {
    dof_handler.clear();
}

template <int dim, int spacedim>
void elastic_rod<dim, spacedim>::third_grid(
        const double &len,
        const double &radius,
        const int &n_of_ref) {

    std::cout << "    length of the rod: " << len << std::endl
            << "    radius of the rod: " << radius << std::endl
            << std::endl << std::endl;

    GridGenerator::hyper_cube(triangulation);
    // 	GridGenerator::cylinder ( triangulation, radius, len );
    // 	GridTools::rotate ( -M_PI / 2 , 1, triangulation );
    // 	const CylindricalManifold<3,3> cyl_manifold (
    // 		2,
    // 		1e-10
    // 		);
    // 
    //     triangulation.set_manifold (0, cyl_manifold);
    //     triangulation.set_all_manifold_ids(0);
    triangulation.refine_global(n_of_ref);
    //     triangulation.set_manifold (0);

    std::cout << "   Number of active cells:       "
            << triangulation.n_active_cells()
            << std::endl << std::endl;


}

template <int dim, int spacedim>
void elastic_rod<dim, spacedim>::setup_system(const bool first_step) {

    if (first_step) {
        dof_handler.distribute_dofs(fe);
        solution_old.reinit(dof_handler.n_dofs());
    }

    hanging_node_constraints.clear();
    DoFTools::make_hanging_node_constraints(
            dof_handler,
            hanging_node_constraints);
    hanging_node_constraints.close();

    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp);
    hanging_node_constraints.condense(dsp);
    sparsity_pattern.copy_from(dsp);
    system_matrix.reinit(sparsity_pattern);

}

template <int dim, int spacedim>
void elastic_rod<dim, spacedim>::assemble_system() {

    
    // 	const FEValuesExtractors::Vector centerline (0);
    // 	const FEValuesExtractors::Vector rotations (spacedim);
    const FEValuesExtractors::Scalar r1(0);
    const FEValuesExtractors::Scalar r2(1);
    const FEValuesExtractors::Scalar r3(2);
    const FEValuesExtractors::Scalar theta1(3);
    const FEValuesExtractors::Scalar theta2(4);
    const FEValuesExtractors::Scalar theta3(5);


    QGauss<dim> quadrature_formula( quadrature );
    FEValues<dim, spacedim> fe_values(fe, quadrature_formula,
            update_values |
            update_gradients |
            update_quadrature_points |
            update_JxW_values);

    QGauss < dim - 1 > face_quadrature_formula( face_quadrature );
    FEFaceValues<dim> fe_face_values(fe, face_quadrature_formula,
            update_values |
            update_gradients |
            update_quadrature_points |
            // 		update_normal_vectors |
            update_JxW_values
            );


    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points = quadrature_formula.size();
    const unsigned int n_face_q_points = face_quadrature_formula.size();


    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double> cell_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);


    material_law.set_old_solution( solution_old );

    typename DoFHandler<dim, spacedim>::active_cell_iterator
    cell = dof_handler.begin_active(),
            endc = dof_handler.end();
    for (; cell != endc; ++cell) // ELEMENT
    {
        cell_matrix = 0;
        cell_rhs = 0;
        fe_values.reinit(cell);


        // GAx, GAy, EA, EIx, EIy, GJ
        material_law.set_coefficients();
        material_law.set_solution(fe_values, n_q_points);


        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) // QUADRATUR 
        {
            const double s_qpoint = fe_values.quadrature_point(q_point)[0];
            material_law.get_solution(q_point, s_qpoint);
            Tensor < 1, 3 > temp;
            temp = cross_product_3d(material_law.get_drds(), material_law.get_n());
            //            material_law.print_member();

            for (unsigned int i = 0; i < dofs_per_cell; ++i) // DOFS(i)
            {
                Tensor < 1, 3 > phi_ir_value;
                Tensor < 1, 3 > phi_ir_grad;
                Tensor < 1, 3 > psi_iR_value;
                Tensor < 1, 3 > psi_iR_grad;

                phi_ir_value[0] = fe_values[r1].value(i, q_point);
                phi_ir_value[1] = fe_values[r2].value(i, q_point);
                phi_ir_value[2] = fe_values[r3].value(i, q_point);

                phi_ir_grad[0] = fe_values[r1].gradient(i, q_point)[0];
                phi_ir_grad[1] = fe_values[r2].gradient(i, q_point)[0];
                phi_ir_grad[2] = fe_values[r3].gradient(i, q_point)[0];

                psi_iR_value[0] = fe_values[theta1].value(i, q_point);
                psi_iR_value[1] = fe_values[theta2].value(i, q_point);
                psi_iR_value[2] = fe_values[theta3].value(i, q_point);

                psi_iR_grad[0] = fe_values[theta1].gradient(i, q_point)[0];
                psi_iR_grad[1] = fe_values[theta2].gradient(i, q_point)[0];
                psi_iR_grad[2] = fe_values[theta3].gradient(i, q_point)[0];


                // RHS: Force and moment densities
                cell_rhs(i) += ((
                        phi_ir_value * loads.f_n_hat(s_qpoint, material_law.get_R())
                        + psi_iR_value * loads.f_m_hat(s_qpoint, material_law.get_R())
                        ) * fe_values.JxW(q_point));

                // term of linearization
                cell_rhs(i) += ((
                        -phi_ir_grad * material_law.get_n()
                        - psi_iR_grad * material_law.get_m()
                        + psi_iR_value * temp
                        ) * fe_values.JxW(q_point));
                // --------------

                for (unsigned int j = 0; j < dofs_per_cell; ++j) // DOFS(j)
                {
                    Tensor < 1, 3 > phi_jr_grad;
                    Tensor < 1, 3 > psi_jR_value;
                    Tensor < 1, 3 > psi_jR_grad;

                    phi_jr_grad[0] = fe_values[r1].gradient(j, q_point)[0];
                    phi_jr_grad[1] = fe_values[r2].gradient(j, q_point)[0];
                    phi_jr_grad[2] = fe_values[r3].gradient(j, q_point)[0];

                    psi_jR_value[0] = fe_values[theta1].value(j, q_point);
                    psi_jR_value[1] = fe_values[theta2].value(j, q_point);
                    psi_jR_value[2] = fe_values[theta3].value(j, q_point);

                    psi_jR_grad[0] = fe_values[theta1].gradient(j, q_point)[0];
                    psi_jR_grad[1] = fe_values[theta2].gradient(j, q_point)[0];
                    psi_jR_grad[2] = fe_values[theta3].gradient(j, q_point)[0];


                    // call side calculation
                    Tensor < 1, 3 > result_1;
                    Tensor < 1, 3 > result_2;
                    Tensor < 1, 3 > result_3;

                    side_calculation(
                            psi_jR_value,
                            psi_jR_grad,
                            phi_jr_grad,
                            result_1,
                            result_2,
                            result_3);

                    // add result_1
                    cell_matrix(i, j) += (
                            (phi_ir_grad * result_1)
                            * fe_values.JxW(q_point));

                    // add result_2
                    cell_matrix(i, j) += (
                            (psi_iR_value * result_2)
                            * fe_values.JxW(q_point));

                    // add result_3
                    cell_matrix(i, j) += (
                            (psi_iR_grad * result_3)
                            * fe_values.JxW(q_point));


                } // end j
            } // end i 
        } // end q


        // RHS: Neumann
        for (unsigned int face_number = 0; face_number < GeometryInfo<dim>::faces_per_cell; ++face_number)
            if (cell->face(face_number)->at_boundary()) // && (cell->face(face_number)->boundary_id() == my_neumann_id)
            {
                /*
                 * cell->face(i)->at_boundary() is always 1 for faces at the boundary
                 * face_number is here 0 (most left) or 1 (most right)
                 */
                fe_face_values.reinit(cell, face_number);
                material_law.set_coefficients();
                material_law.set_solution(fe_face_values, n_face_q_points);

                for (unsigned int q_point = 0; q_point < n_face_q_points; ++q_point) // here mostly one
                {
                    const double face_s_qpoint = fe_face_values.quadrature_point(q_point)[0]; // s
                    material_law.get_solution(q_point, face_s_qpoint);

                    for (unsigned int i = 0; i < dofs_per_cell; ++i) {

                        Tensor < 1, 3 > face_phi_ir_value;
                        Tensor < 1, 3 > face_psi_iR_value;

                        face_phi_ir_value[0] = fe_face_values[r1].value(i, q_point);
                        face_phi_ir_value[1] = fe_face_values[r2].value(i, q_point);
                        face_phi_ir_value[2] = fe_face_values[r3].value(i, q_point);

                        face_psi_iR_value[0] = fe_face_values[theta1].value(i, q_point);
                        face_psi_iR_value[1] = fe_face_values[theta2].value(i, q_point);
                        face_psi_iR_value[2] = fe_face_values[theta3].value(i, q_point);

                        int entry = 99;
                        for (int i = 0; i < 3; i++) {
                            if (face_phi_ir_value[i] == 1) {
                                entry = i;
                            }
                            if (face_psi_iR_value[i] == 1) {
                                entry = i + 3;
                            }
                        }
                        if (loads.ask_for_boundary(face_number, entry) == true) {

                            /*
                             * adding prescribed force terms 
                             */
                            cell_rhs(i) += (
                                    face_phi_ir_value * loads.f_n_pre(material_law.get_R())
                                    + face_psi_iR_value * loads.f_m_pre(material_law.get_R())
                                    ) * fe_face_values.JxW(q_point);

                        } else {

                            // influences only for dirichlet oundary != 0 
                            cell_rhs(i) += ((
                                    face_phi_ir_value * material_law.get_n()
                                    + face_psi_iR_value * material_law.get_m()
                                    ) * fe_face_values.JxW(q_point));

                            
                            /* 
                             * adding unknowns to linear system 
                             */
                            for (unsigned int j = 0; j < dofs_per_cell; ++j) {

                                Tensor < 1, 3 > face_phi_jr_grad;
                                Tensor < 1, 3 > face_psi_jR_value;
                                Tensor < 1, 3 > face_psi_jR_grad;

                                face_phi_jr_grad[0] = fe_face_values[r1].gradient(j, q_point)[0];
                                face_phi_jr_grad[1] = fe_face_values[r2].gradient(j, q_point)[0];
                                face_phi_jr_grad[2] = fe_face_values[r3].gradient(j, q_point)[0];

                                face_psi_jR_value[0] = fe_face_values[theta1].value(j, q_point);
                                face_psi_jR_value[1] = fe_face_values[theta2].value(j, q_point);
                                face_psi_jR_value[2] = fe_face_values[theta3].value(j, q_point);

                                face_psi_jR_grad[0] = fe_face_values[theta1].gradient(j, q_point)[0];
                                face_psi_jR_grad[1] = fe_face_values[theta2].gradient(j, q_point)[0];
                                face_psi_jR_grad[2] = fe_face_values[theta3].gradient(j, q_point)[0];


                                Tensor < 1, 3 > result_1;
                                Tensor < 1, 3 > result_2;
                                Tensor < 1, 3 > result_3;

                                side_calculation(
                                        face_psi_jR_value,
                                        face_psi_jR_grad,
                                        face_phi_jr_grad,
                                        result_1,
                                        result_2,
                                        result_3);

                                cell_matrix(i, j) += (
                                        (face_phi_ir_value * result_1)
                                        * fe_face_values.JxW(q_point));


                                cell_matrix(i, j) += (
                                        (face_psi_iR_value * result_3)
                                        * fe_face_values.JxW(q_point));

                            } // j
                        } // if: rhs known - rhs unknown
                    } // i  
                } // q_point  
            } // for faces - if boundary 



        cell->get_dof_indices(local_dof_indices);
        for (unsigned int i = 0; i < dofs_per_cell; ++i) {
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
                system_matrix.add(local_dof_indices[i],
                    local_dof_indices[j],
                    cell_matrix(i, j));

            system_rhs(local_dof_indices[i]) += cell_rhs(i);
        }
    } // end element 


    hanging_node_constraints.condense(system_matrix);
    hanging_node_constraints.condense(system_rhs);

    // Idee: alle Masken vorbereiten und Boundary IDs indivduell setzen
    // 12 Optionen (links, rechts, 6 Komponenten)
    std::map<types::global_dof_index, double> boundary_values;
    FEValuesExtractors::Scalar used_mask[] = {r1, r2, r3, theta1, theta2, theta3};
    for (int i = 0; i < 6; i++) {
        // most left boundary
        
        if (loads.ask_for_boundary(0, i) == false)
            VectorTools::interpolate_boundary_values(
                dof_handler, 0,
//                ZeroFunction<spacedim > (6),
                ConstantFunction<spacedim > ( loads.get_dirichlet_left(i) , 6),
                boundary_values,
                fe.component_mask(used_mask[i])
                );

        // most right boundary
        if (loads.ask_for_boundary(1, i) == false)
            VectorTools::interpolate_boundary_values(
                dof_handler, 1,
//                ZeroFunction<spacedim > (6),
                ConstantFunction<spacedim > (loads.get_dirichlet_right(i) , 6),
                boundary_values,
                fe.component_mask(used_mask[i])
                );
    }

    MatrixTools::apply_boundary_values(
            boundary_values, system_matrix, solution, system_rhs);
    
    
    return;
}

template <int dim, int spacedim>
void elastic_rod<dim, spacedim>::solve(bool first_call ) {
    SolverControl solver_control(10000, 1e-12);
    SolverCG<> cg(solver_control);


//    std::cout << " matrix to file: ";
//    std::ofstream of("sys_mat.txt");
//    system_matrix.print_formatted(of);
//    std::cout << " - done - ";
    std::cout << "Solve: ";

    try {
        // Direct
        Timer timer;
        timer.start();
        SparseDirectUMFPACK A_direct;
        A_direct.initialize(system_matrix);
        A_direct.vmult(solution, system_rhs);
        timer.stop();
        std::cout << "direct - done (" << timer() << "s)"
                << std::endl;
    } catch (...) {
        std::cerr
                << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
        std::cerr
                << "Could not use direct solver!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

        // CG
        Timer timer_cg;
        timer_cg.start();
        //     PreconditionSSOR<> preconditioner;
        //     preconditioner.initialize(system_matrix, 1.2);
        // 
        cg.solve(
                system_matrix,
                solution,
                system_rhs,
                PreconditionIdentity());
        timer_cg.stop();
        std::cout << "cg - done (" << timer_cg() << "s)"
                << std::endl;
    }

    hanging_node_constraints.distribute(solution);


    // evaluation
    residual.push_back(system_rhs.l2_norm());

    
    // update solution
    double step_length = get_newton_step_length( first_call );
    std::cout <<"  alpha: "<<  step_length << " - "; 
    alphas.push_back( step_length );
    solution_old.add(step_length, solution);
}

/**
 * Output:
 * 1) only components
 * 2) vector values (does not work)
 */
template <int dim, int spacedim>
void elastic_rod<dim, spacedim>::output_results(const unsigned int cycle) const {

    std::ostringstream oss;
    oss << "sol/solution-" << cycle << ".vtk";

    std::string filename = oss.str();
    //    if (cycle < 10) {
    //        filename += ('0' + cycle);
    //    } else {
    //        filename += ( '1' + cycle%10 );
    //    }        
    //    Assert(cycle < 100, ExcInternalError());
    //
    //    filename += ".vtk";
    std::ofstream output(filename.c_str());

    std::vector<std::string> solution_names;
    solution_names.push_back("rx");
    solution_names.push_back("ry");
    solution_names.push_back("rz");
    solution_names.push_back("theta1");
    solution_names.push_back("theta2");
    solution_names.push_back("theta3");

    DataOut<dim, DoFHandler<dim, spacedim >> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution_old, solution_names);
    data_out.build_patches();
    data_out.write_vtk(output);


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
double elastic_rod<dim, spacedim>::get_newton_step_length(bool first_call) {
    /* find optimal step length: trust region or back tracking line search */

    /*
     *  backtracking:
     * - choose a=1
     * - if r(a=1) << r(a=0) then choose a=1
     * - else check r(a=2/3) << r(a=0) if so then choos
     * - else check r(a=(2/3)^2 ) << r( ... ) do so ..
     * - ...
     * 
     * 2/3 can be chosen different - look for Wolfe and Armijo-Goldstein condition
     * 
     */
    
    double alpha = 1.0;
    
    if ( global_alpha_off ) {
        res_check_all.push_back( compute_residual( alpha ) );
        return alpha;
    }
    
//    if ( first_call ) {
//        res_check_all.push_back( compute_residual( alpha ) );
//        return alpha;
//    }
    
    double res_zero;
    res_zero = residual.back();
    
//    if (first_call) { 
//        res_zero = residual.back();
////        return alpha;
//    } else {
//        res_zero = res_check_all.back();
//    }
    
    double res_check;
    while ( alpha > 0.1 ) {
    
        
        /* compute residual for alpha */
        res_check = compute_residual( alpha );

        /* check residual */
        if ( res_check * 1.0 < res_zero ) {
            /* good enough */
            res_check_all.push_back( res_check );
            return alpha;
        }

        /* update alpha */
        alpha *= 2.0 / 3.0;        
    }
    
    res_check_all.push_back( res_check );
    return alpha;
}

template <int dim, int spacedim>
void elastic_rod<dim, spacedim>::side_calculation(
        const Tensor < 1, 3 > &psi_jR_value,
        const Tensor < 1, 3 > &psi_jR_grad,
        const Tensor < 1, 3 > &phi_jr_grad,
        Tensor < 1, 3 > &result_1,
        Tensor < 1, 3 > &result_2,
        Tensor < 1, 3 > &result_3
        ) {

    // forces and moments
    const Tensor < 1, 3 > n_vec = material_law.get_n();
    const Tensor < 1, 3 > m_vec = material_law.get_m();
    const Tensor < 1, 3 > dr_ds_old = material_law.get_drds();


    // helpers
    Tensor < 1, 3 > temp_a;
    Tensor < 1, 3 > temp_b;
    Tensor < 1, 3 > temp_c;

    temp_a = cross_product_3d(dr_ds_old, psi_jR_value);
    temp_a += phi_jr_grad;

    // 1st part:
    temp_b = cross_product_3d(n_vec, psi_jR_value);
    temp_b *= -1;
    temp_c = material_law.get_R_Cmn_RT(1, 1) * temp_a;
    temp_b += temp_c;
    temp_c = material_law.get_R_Cmn_RT(1, 2) * psi_jR_grad;
    temp_b += temp_c;
    result_1 = temp_b;

    // 2nd part
    temp_c = cross_product_3d(dr_ds_old, temp_b);
    temp_b = cross_product_3d(n_vec, phi_jr_grad);
    result_2 = temp_b - temp_c;

    // 3rd part 
    temp_b = cross_product_3d(m_vec, psi_jR_value);
    temp_b *= -1;
    temp_c = material_law.get_R_Cmn_RT(2, 1) * temp_a;
    temp_b += temp_c;
    temp_c = material_law.get_R_Cmn_RT(2, 2) * psi_jR_grad;
    temp_b += temp_c;
    result_3 = temp_b;

}

template <int dim, int spacedim>
double elastic_rod<dim, spacedim>::compute_residual(double alpha) {

    Vector<double> residual;
    residual.reinit(dof_handler.n_dofs());
    Vector<double> evaluation_point = solution_old;
    evaluation_point.add(alpha, solution);

    // 	const FEValuesExtractors::Vector centerline (0);
    // 	const FEValuesExtractors::Vector rotations (spacedim);
    const FEValuesExtractors::Scalar r1(0);
    const FEValuesExtractors::Scalar r2(1);
    const FEValuesExtractors::Scalar r3(2);
    const FEValuesExtractors::Scalar theta1(3);
    const FEValuesExtractors::Scalar theta2(4);
    const FEValuesExtractors::Scalar theta3(5);


    QGauss<dim> quadrature_formula( quadrature );
    FEValues<dim, spacedim> fe_values(fe, quadrature_formula,
            update_values |
            update_gradients |
            update_quadrature_points |
            update_JxW_values);

    QGauss < dim - 1 > face_quadrature_formula( face_quadrature );
    FEFaceValues<dim> fe_face_values(fe, face_quadrature_formula,
            update_values |
            update_gradients |
            update_quadrature_points |
            // 		update_normal_vectors |
            update_JxW_values
            );


    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points = quadrature_formula.size();
    const unsigned int n_face_q_points = face_quadrature_formula.size();


    Vector<double> cell_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    material_law.set_old_solution(evaluation_point);

    typename DoFHandler<dim, spacedim>::active_cell_iterator
    cell = dof_handler.begin_active(),
            endc = dof_handler.end();
    for (; cell != endc; ++cell) // ELEMENT
    {
        cell_rhs = 0;
        fe_values.reinit(cell);


        // GAx, GAy, EA, EIx, EIy, GJ
        material_law.set_coefficients();
        material_law.set_solution(fe_values, n_q_points);


        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) // QUADRATUR 
        {
            const double s_qpoint = fe_values.quadrature_point(q_point)[0];
            material_law.get_solution(q_point, s_qpoint);
            Tensor < 1, 3 > temp;
            temp = cross_product_3d(material_law.get_drds(), material_law.get_n());
            //            material_law.print_member();

            for (unsigned int i = 0; i < dofs_per_cell; ++i) // DOFS(i)
            {
                Tensor < 1, 3 > phi_ir_value;
                Tensor < 1, 3 > phi_ir_grad;
                Tensor < 1, 3 > psi_iR_value;
                Tensor < 1, 3 > psi_iR_grad;

                phi_ir_value[0] = fe_values[r1].value(i, q_point);
                phi_ir_value[1] = fe_values[r2].value(i, q_point);
                phi_ir_value[2] = fe_values[r3].value(i, q_point);

                phi_ir_grad[0] = fe_values[r1].gradient(i, q_point)[0];
                phi_ir_grad[1] = fe_values[r2].gradient(i, q_point)[0];
                phi_ir_grad[2] = fe_values[r3].gradient(i, q_point)[0];

                psi_iR_value[0] = fe_values[theta1].value(i, q_point);
                psi_iR_value[1] = fe_values[theta2].value(i, q_point);
                psi_iR_value[2] = fe_values[theta3].value(i, q_point);

                psi_iR_grad[0] = fe_values[theta1].gradient(i, q_point)[0];
                psi_iR_grad[1] = fe_values[theta2].gradient(i, q_point)[0];
                psi_iR_grad[2] = fe_values[theta3].gradient(i, q_point)[0];


                // RHS: Force and moment densities
                cell_rhs(i) += ((
                        phi_ir_value * loads.f_n_hat(s_qpoint, material_law.get_R())
                        + psi_iR_value * loads.f_m_hat(s_qpoint, material_law.get_R())
                        ) * fe_values.JxW(q_point));

                // term of linearization
                cell_rhs(i) += ((
                        -phi_ir_grad * material_law.get_n()
                        - psi_iR_grad * material_law.get_m()
                        + psi_iR_value * temp
                        ) * fe_values.JxW(q_point));
            } // end i 
        } // end q


        // RHS: Neumann
        for (unsigned int face_number = 0; face_number < GeometryInfo<dim>::faces_per_cell; ++face_number)
            if (cell->face(face_number)->at_boundary()) // && (cell->face(face_number)->boundary_id() == my_neumann_id)
            {
                /*
                 * cell->face(i)->at_boundary() is always 1 for faces at the boundary
                 * face_number is here 0 (most left) or 1 (most right)
                 */
                fe_face_values.reinit(cell, face_number);
                material_law.set_coefficients();
                material_law.set_solution(fe_face_values, n_face_q_points);

                for (unsigned int q_point = 0; q_point < n_face_q_points; ++q_point) // here mostly one
                {
                    const double face_s_qpoint = fe_face_values.quadrature_point(q_point)[0]; // s
                    material_law.get_solution(q_point, face_s_qpoint);

                    for (unsigned int i = 0; i < dofs_per_cell; ++i) {

                        Tensor < 1, 3 > face_phi_ir_value;
                        Tensor < 1, 3 > face_psi_iR_value;

                        face_phi_ir_value[0] = fe_face_values[r1].value(i, q_point);
                        face_phi_ir_value[1] = fe_face_values[r2].value(i, q_point);
                        face_phi_ir_value[2] = fe_face_values[r3].value(i, q_point);

                        face_psi_iR_value[0] = fe_face_values[theta1].value(i, q_point);
                        face_psi_iR_value[1] = fe_face_values[theta2].value(i, q_point);
                        face_psi_iR_value[2] = fe_face_values[theta3].value(i, q_point);

                        int entry = 99;
                        for (int i = 0; i < 3; i++) {
                            if (face_phi_ir_value[i] == 1) {
                                entry = i;
                            }
                            if (face_psi_iR_value[i] == 1) {
                                entry = i + 3;
                            }
                        }
                        if (loads.ask_for_boundary(face_number, entry) == true) {

                            /*
                             * adding prescribed force terms 
                             */
                            cell_rhs(i) += (
                                    face_phi_ir_value * loads.f_n_pre(material_law.get_R())
                                    + face_psi_iR_value * loads.f_m_pre(material_law.get_R())
                                    ) * fe_face_values.JxW(q_point);

                        } else {

                            // influences only for dirichlet oundary != 0 
                            cell_rhs(i) += ((
                                    face_phi_ir_value * material_law.get_n()
                                    + face_psi_iR_value * material_law.get_m()
                                    ) * fe_face_values.JxW(q_point));
                        } // if: rhs known - rhs unknown
                    } // i  
                } // q_point  
            } // for faces - if boundary 

        cell->get_dof_indices(local_dof_indices);
        for (unsigned int i = 0; i < dofs_per_cell; ++i) {
            residual(local_dof_indices[i]) += cell_rhs(i);
        }
    } // end element 

    hanging_node_constraints.condense(system_matrix);
    hanging_node_constraints.condense(residual);

    // Idee: alle Masken vorbereiten und Boundary IDs indivduell setzen
    // 12 Optionen (links, rechts, 6 Komponenten)
    std::map<types::global_dof_index, double> boundary_values;
    FEValuesExtractors::Scalar used_mask[] = {r1, r2, r3, theta1, theta2, theta3};
    for (int i = 0; i < 6; i++) {
        // most left boundary
        if (loads.ask_for_boundary(0, i) == false)
            VectorTools::interpolate_boundary_values(
                dof_handler, 0,
                ZeroFunction<spacedim > (6),
                boundary_values,
                fe.component_mask(used_mask[i])
                );

        // most right boundary
        if (loads.ask_for_boundary(1, i) == false)
            VectorTools::interpolate_boundary_values(
                dof_handler, 1,
                ZeroFunction<spacedim > (6),
                boundary_values,
                fe.component_mask(used_mask[i])
                );
    }

    MatrixTools::apply_boundary_values(
            boundary_values, system_matrix, solution, residual);

    return residual.l2_norm();
}

template <int dim, int spacedim>
void elastic_rod<dim, spacedim>::run() {

    bool first_call = true;
    residual.resize(0);
    double old_residual = 0;

    /* geometry, mesh, triangulation, grid */
    third_grid(1, 0.2, global_refinement);
    
    /* ... */
    const double my_error = 1.0e-9;
    
    
    // newton itarationen
    for (unsigned int newton_count = 0; newton_count < max_newton_iter ; ++newton_count) {
        std::cout << "newton: " << newton_count;

        /* initialize linear system */
        setup_system(first_call);

        std::cout << " - DOFS: "
                << dof_handler.n_dofs() << " - ";

        /* assemble linear system */
        assemble_system();

        /* solve linear system */
        solve( first_call );

        /* write solution to file */
        output_results(newton_count);

        /* check newton accuracy and progress*/
//        if (std::abs(old_residual / residual.back() - 1) < 0.01) {
//            std::cout << " ---  enough ---" << '\n' << '\n';
//            break;
//        }
        
        if ( !first_call && residual.back() < 100.0 )
            loads.increase_load();
        
        if ( !first_call && residual.back() < my_error ) {
            if ( loads.increase_load() ) {}
            else {
                std::cout       << "\n\n >>> good enough err:  " 
                                << residual.back() <<  " <<< \n\n";
                break;
            }
        }
        
        const double my_improvement = abs(residual[residual.size()-2] - residual.back());
        if ( !first_call && my_improvement < my_error  ) {
                if ( loads.increase_load() ) {}
                else {
                std::cout       << "\n\n >>> no improvement any more:  "
                                << residual.back() <<  " <<< \n\n";
                break;
                }
         }
        
        if (!first_call && residual.back() > old_residual) {
            std::cout << " -!- residual grows -!- ";
        }

        /* update old_residual */
        old_residual = residual.back();
        first_call = false;
    }

    std::cout
            << std::endl << " ---------- " << std::endl
            << "  num of newton iterations: " << residual.size() << std::endl;

//    std::vector<double>::iterator item = residual.begin();
//    for (; item != residual.end(); ++item) {
//        std::cout << "     " << *item << std::endl;
//    }
    for (int i=0; i< residual.size() ;i++)
    {
        std::cout 
                << "\t" << "step-" << i << " ------- " << '\n' 
                << "\t\t" << "residual (x+dx): " << res_check_all[i] << '\n' 
//                << "\t\t" << "residual: " << residual[i] << '\n' 
                << "\t\t" << "alpha: " << alphas[i] << std::endl;
    }
    std::cout << std::endl << " ---------- " << std::endl;

}