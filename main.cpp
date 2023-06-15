#include <functional>
#include <stdexcept>
#include <iostream>
#include <numeric>
#include <sstream>
#include <fstream>
#include <vector>
#include <memory>
#include <cmath>
#include <regex>


#include "Complex.h"
#include "Vector.h"
#include "BlockMatrix.h"


#define MKL_Complex16 Co <double>


#include <mkl.h>
#include <omp.h>


#include "Specializations.h"
#include "Quadrature.h"
#include "Mesh.h"
//#include "Solver.h"


using Scalar = double;
using Length = MKL_INT;


double arr [6][6] =
{
    { 3.966261, 2.133467, 1.216731, 5.262742, 4.057779, 8.091164, },
    { 2.133467, 5.288848, 6.238135, 2.130408, 6.813768, 8.245826, },
    { 1.216731, 6.238135, 5.880156, 0.952626, 3.437681, 6.434359, },
    { 5.262742, 2.130408, 0.952626, 6.579653, 2.719752, 2.975403, },
    { 4.057779, 6.813768, 3.437681, 2.719752, 4.190280, 5.671886, },
    { 8.091164, 8.245826, 6.434359, 2.975403, 5.671886, 7.223314, },
};


int main ()
{
    BlockSystem <Scalar, Length> block_system("system");

    block_system.set_block_size(1);

    block_system.set_matrix_size(6);

    block_system.set_column_size(2, 6);

    block_system.prepare_system([&] (Length col, Length row, Scalar & val)
    {
        val = arr[row][col];
    });

    std::cout << block_system.matrix_to_string() << std::endl;

    block_system.prepare_column([&] (Length col, Length row, Scalar & val)
    {
        val = 1.0 + (Scalar) row;
    });

    std::cout << block_system.column_to_string() << std::endl;

    block_system.prepare_solution();

    std::cout << block_system.column_to_string() << std::endl;

    return 0;
}
