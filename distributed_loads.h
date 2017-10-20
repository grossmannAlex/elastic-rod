/* 
 * File:   distributed_loads.h
 * Author: alex
 *
 * Created on October 13, 2017, 2:20 PM
 */

#ifndef DISTRIBUTED_LOADS_H
#define	DISTRIBUTED_LOADS_H


#include <deal.II/lac/vector.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/diagonal_matrix.h>
#include <deal.II/base/exceptions.h>


using namespace dealii;


class distributed_loads {
public:
    distributed_loads();
    distributed_loads(const distributed_loads& orig);
    virtual ~distributed_loads();
    
    // distributed loads
    Tensor<1, 3> f_n_hat( const double &s, const Tensor < 2, 3 > R );
    Tensor<1, 3> f_m_hat( const double &s, const Tensor < 2, 3 > R );
    
    // neumann
    Tensor<1, 3> f_n_pre( const Tensor < 2, 3 > R);
    Tensor<1, 3> f_m_pre( const Tensor < 2, 3 > R);
    
    double get_dirichlet_left( int entry_num );
    double get_dirichlet_right( int entry_num );
    
    bool ask_for_boundary(
        const unsigned int i,
        const unsigned int j);
private:

};

#endif	/* DISTRIBUTED_LOADS_H */
