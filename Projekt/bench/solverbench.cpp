//#define NODEBUG

#include "mpihandler.hpp"
#include "distellpackmatrix.hpp"
#include "solverbench.hpp"

int main()
{
    auto mat = Icarus::construct_model_matrix<Icarus::DistEllpackMatrix<double>>(2);
    Icarus::print_sliced_object(mat);
    return 0;
}
