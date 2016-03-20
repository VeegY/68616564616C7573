#define NODEBUG
#define LOGGING_LEVEL 2

#include "mpihandler.hpp"
#include "distellpackmatrix.hpp"
#include "solverbench.hpp"

int main()
{
    //check model matrix
    //auto mat = Icarus::construct_model_matrix<Icarus::DistEllpackMatrixGpu<double>>(4);
    //Icarus::print_sliced_object(mat);
    Icarus::SolverBench<Icarus::DistEllpackMatrix<double>> benchmark({1},{128});
    benchmark.run();
    benchmark.print_results(std::cout);
    return 0;
}
