/* ---------------------------------------------------------------------
 *
 * Copyright (C) 
 *
 * Author: Alexander Grossmann, University of Erlangen-Nuernberg, 2017
 */



#include <fstream>
#include <iostream>
#include <stdlib.h>     /* abs */
#include <math.h>


// classes
#include "elastic_rod.cc"




int main() {
    using namespace dealii;
    
    try {


//        StepRod::
        elastic_rod < 1, 1 > elastic_problem_3d;
        elastic_problem_3d.run();
        
    } catch (std::exception &exc) {
        std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
        std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

        return 1;
    } catch (...) {
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
