/* 
 * File:   elastic_rod.h
 * Author: alex
 *
 * Created on October 13, 2017, 2:24 PM
 */

#ifndef ELASTIC_ROD_H
#define	ELASTIC_ROD_H

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
// #include <deal.II/lac/constraint_matrix.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/base/timer.h>


#include <fstream>
#include <iostream>
#include <sstream>
#include <stdlib.h>     /* abs */
#include <math.h>


#include "distributed_loads.cc"
#include "constitutive_law.cc"

    
using namespace dealii;

template <int dim, int spacedim>
class elastic_rod {
public:
    elastic_rod();
    ~elastic_rod();
    void run();

private:
    void third_grid(
            const double &len,
            const double &radius,
            const int &n_of_ref);
    void setup_system(const bool first_step);
    void assemble_system( );
    void solve( bool first_step );
    void refine_grid();
    void output_results(const unsigned int cycle) const;


    void side_calculation(
            const Tensor < 1, 3 > &psi_jR_value,
            const Tensor < 1, 3 > &psi_jR_grad,
            const Tensor < 1, 3 > &phi_jr_grad,
            Tensor < 1, 3 > &result_1,
            Tensor < 1, 3 > &result_2,
            Tensor < 1, 3 > &result_3
            );


    /* newton iterations methods */
    double get_newton_step_length( bool first_call );
    double compute_residual (double alpha);
    
    
    constitutive_law material_law;
    distributed_loads loads;

    Triangulation<dim, spacedim> triangulation;
    DoFHandler<dim, spacedim> dof_handler;
    FESystem<dim, spacedim> fe;

    ConstraintMatrix hanging_node_constraints;
    SparsityPattern sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double> solution;
    Vector<double> solution_old;
    Vector<double> system_rhs;
    
    std::vector<double> residual;
    std::vector<double> res_check_all;
    std::vector<double> alphas;

};


#endif	/* ELASTIC_ROD_H */

