#define NODEBUG

#include "mpihandler.hpp"
#include "distellpackmatrix.hpp"
#include "solverbench.hpp"

int main()
{
    //check model matrix
    //auto mat = Icarus::construct_model_matrix<Icarus::DistEllpackMatrix<double>>(3);
    //Icarus::print_sliced_object(mat);
    
    Icarus::SolverBench<Icarus::DistEllpackMatrix<double>> benchmark({1,2,4,8,16},{8,16,32,64,128,256,512,1024,2048,4096});
    benchmark.run();
    benchmark.print_results(std::cout);
    return 0;
}
