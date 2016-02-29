#include "../src/include/slicedvector.hpp"
#include "../src/include/distellpackmatrix.hpp"
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <limits>
/*
*test fuer Aequilibrierungs-Vorkonditionierer
*
*/
int equipc_test()
{
    const size_t N = 1000;
    const size_t maxrow = 10;
    //create random matrix
    srand(static_cast <unsigned> (time(0)));
    Icarus::DistEllpackMatrix<double> mat1(N), equi(N);

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
        for (size_t j(0); j <= rowlen; j++)
        {
            val = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);  //all values positive
            colind = rand() % N;
            mat1.sequential_fill(colind, val);
        }
        mat1.end_of_row();
    }
    if (!mat1.is_filled()) LOG_ERROR("filling of random Matrix failed");
    equi = mat1.precond_equi();

    //initialize testvector
    Icarus::SlicedVector<double> rhs(N);
    for (size_t i(fron); i <= lron; i++)
    {
        rhs.set_global(i, 1.0);
    }
    mat1.mult_vec(rhs, rhs);
    equi.mult_vec(rhs, rhs);
    double checktol;
    checktol = std::numeric_limits<double>::epsilon() * 10.0;
    for (size_t i(fron); i < lron; i++)
    {
        if (std::abs(rhs.get_global(i) - 1.0) > checktol)
            LOG_ERROR("equi pc failed; difference to 1.0: ", std::abs(rhs.get_global(i) - 1.0));
    }
    return 0;
}

int main()
{
    return equipc_test();
}
