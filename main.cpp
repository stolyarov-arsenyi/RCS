#include <functional>
#include <filesystem>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <fstream>
#include <vector>
#include <memory>
#include <cmath>
#include <regex>


#include "Complex.h"
#include "Vector.h"
#include "BlockSystem.h"


#define MKL_Complex16 Co <double>


#include <mkl.h>
#include <omp.h>


#include "Specializations.h"
#include "Solver.h"


using Scalar = double;
using Length = MKL_INT;


int main (int argc, char ** argv)
{
    if (argc != 3)
    {
        std::cout << "Wrong number of arguments" << std::endl;

        return 0;
    }

    std::string conf_name = argv[1];
    std::string mesh_name = argv[2];

    Solver <Scalar, Length> solver(conf_name, mesh_name);

    solver.init();
    solver.solve();

    return 0;
}
