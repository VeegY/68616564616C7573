#include "../src/include/slicedvectorgpu.hpp"
#include "../src/include/bicgstabsolver.hpp"
#include "../src/include/distellpackmatrixgpu.hpp"
#include "../src/include/slicedvector.hpp"
#include "../src/include/distellpackmatrix.hpp"

#include <cmath>

#define LOGGING_LEVEL 3
#include "../src/include/logger.hpp"

int main()
{
    Icarus::DistEllpackMatrixGpu<double> matrix(3);
    Icarus::SlicedVectorGpu<double> rhs(3);
    Icarus::SlicedVectorGpu<double> sol(3);
    Icarus::BiCgStabSolver<Icarus::DistEllpackMatrixGpu<double>> solver(matrix, rhs);
    matrix.prepare_sequential_fill(1);
    matrix.sequential_fill(0, 1.0);
    matrix.end_of_row();
    matrix.sequential_fill(1, 1.0);
    matrix.end_of_row();
    matrix.sequential_fill(2, 1.0);
    matrix.end_of_row();

    rhs.set_global(0, 1.0);
    rhs.set_global(1, 1.0);
    rhs.set_global(2, 1.0);
    sol.clear();
    sol.set_local(0, 1.0);
    sol.set_local(1, 0.0);
    sol.set_local(2, 0.0);
    solver.solve(sol);

    return 0;
}
