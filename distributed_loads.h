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

using namespace dealii;


class distributed_loads {
public:
    distributed_loads();
    distributed_loads(const distributed_loads& orig);
    virtual ~distributed_loads();
    
    Tensor<1, 3> f_n_hat( const double &s );
    Tensor<1, 3> f_m_hat( const double &s );
private:

};

#endif	/* DISTRIBUTED_LOADS_H */

