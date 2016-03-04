#include "../src/include/slicedvector.hpp"
#include "../src/include/distellpackmatrix.hpp"
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <limits>
/*
*test fuer jacobi-vorkonditionierer
*
*/
int jacpc_test(size_t N, size_t maxrow)
{
    //create random matrix
    srand(static_cast <unsigned> (time(0)));
    Icarus::DistEllpackMatrix<double> mat1(N), jac(N);
    Icarus::SlicedVector<double> diag_entries(N), rhs(N), res(N);


    //fill mat1 randomly;
    size_t fron, lron;
    fron = mat1.first_row_on_node();
    lron = fron + mat1.get_dim_local() - 1;
    size_t rowlen, colind;
    double val;
    mat1.prepare_sequential_fill(maxrow);
    for (size_t i(fron); i <= lron; i++)
    {
        rowlen = rand() % maxrow;
        val = 2.0 * (static_cast <double> (rand()) / static_cast <double> (RAND_MAX)) - 1.0;  //values between -1 and 1
        diag_entries.set_global(i, val);
        mat1.sequential_fill(i, val);
        for (size_t j(0); j < rowlen; j++)
        {
            val = 2.0 * (static_cast <double> (rand()) / static_cast <double> (RAND_MAX)) - 1.0;  //values between -1 and 1
            colind = rand() % N;
            while (colind==i) colind = rand() % N;
            mat1.sequential_fill(colind, val);
        }
        mat1.end_of_row();
    }
    if (!mat1.is_filled()) LOG_ERROR("filling of random Matrix failed");
    jac = mat1.precond_jacobi();
    for (size_t i(0); i < N; i++)
    {
        val = 2.0 * (static_cast <double> (rand()) / static_cast <double> (RAND_MAX)) - 1.0;  //values between -1 and 1
        rhs.set_global(i, val);
    }
    mat1.mult_vec(rhs, rhs);
    jac.mult_vec(rhs, res);
    double checktol = std::numeric_limits<double>::epsilon()*10;
    for (size_t i(fron); i <= lron; i++)
    {
        if (std::abs((rhs.get_global(i)/diag_entries.get_global(i))-res.get_global(i)) > checktol*std::abs(res.get_global(i)))
            LOG_ERROR("jacobi pc failed;  result: ", res.get_global(i), " ;  error: ",
                       std::abs((rhs.get_global(i)/diag_entries.get_global(i))-res.get_global(i)));
    }
    return 0;
}

int main()
{
    return jacpc_test(100000, 30);
}
