#include <functional>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <memory>
#include <cmath>
#include <regex>


#include "Complex.h"
#include "BlockMatrix.h"


#define MKL_Complex16 Co <double>


#include <mkl.h>
#include <omp.h>


#include "Specializations.h"
#include "Quadrature.h"
#include "Mesh.h"
#include "Solver.h"


using Scalar = MKL_Complex16;
using Length = MKL_INT;


int main ()
{


    BlockMatrix <Scalar, Length> block_matrix(6, 6);

    return 0;
}
