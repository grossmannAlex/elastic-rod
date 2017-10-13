/* 
 * File:   distributed_loads.cc
 * Author: alex
 * 
 * Created on October 13, 2017, 2:20 PM
 */

#include "distributed_loads.h"



distributed_loads::distributed_loads() {
}

distributed_loads::distributed_loads(const distributed_loads& orig) {
}

distributed_loads::~distributed_loads() {
}

Tensor < 1, 3> distributed_loads::f_n_hat(const double &s) {
    using namespace dealii;
    const Tensor<1, 3> result({0, 0, 1});
    return result;
}

Tensor < 1, 3> distributed_loads::f_m_hat(const double &s) {
    const Tensor<1, 3> result({0, 0, 0});
    return result;
}

