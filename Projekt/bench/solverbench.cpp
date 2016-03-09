#define NODEBUG

#include "mpihandler.hpp"
#include "distellpackmatrix.hpp"
#include "solverbench.hpp"

int main()
{
    //auto mat = Icarus::construct_model_matrix<Icarus::DistEllpackMatrix<double>>(3);
    //Icarus::print_sliced_object(mat);
    
    Icarus::SolverBench<Icarus::DistEllpackMatrix<double>> benchmark(1,3,100,150);
    benchmark.run();
    benchmark.print_results(std::cout);
    return 0;
}
