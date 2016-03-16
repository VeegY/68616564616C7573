#define NODEBUG

#include "mpihandler.hpp"
#include "distellpackmatrixgpu.hpp"
#include "solverbench.hpp"

int main()
{
    //check model matrix
    //auto mat = Icarus::construct_model_matrix<Icarus::DistEllpackMatrixGpu<double>>(4);
    //Icarus::print_sliced_object(mat);
    Icarus::SolverBench<Icarus::DistEllpackMatrixGpu<double>> benchmark({1,2,4,8,16},{8,16,32,64,128,256,512,1024, 2048});
    benchmark.run();
    benchmark.print_results(std::cout);
    return 0;
}
