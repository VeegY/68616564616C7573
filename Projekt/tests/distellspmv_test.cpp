#include "../src/include/fullvector.hpp"
#include "../src/include/slicedvector.hpp"
#include "../src/include/distellpackmatrix.hpp"
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <limits>
/*
*test fuer Matrix-Vector Multiplikation im distELLPack Format
*
*/
int distEllSpmv_test(const size_t N, const size_t maxrow)
{
    //create random matrix
    srand(static_cast <unsigned> (time(0)));
    Icarus::DistEllpackMatrix<double> mat1(N);
    Icarus::FullVector<double> res(N);
    res.clear();
    Icarus::SlicedVector<double> rhs(N);
    //fill mat1 randomly;
    size_t fron, lron;
    fron = mat1.first_row_on_node();
    lron = fron + mat1.get_dim_local() - 1;
    size_t rowlen, colind;
    double val;
    mat1.prepare_sequential_fill(maxrow);
    
    for(size_t i=fron; i <= lron; i++)
        rhs.set_global(i, (1./N)*(i+1));

    for (size_t i(fron); i <= lron; i++)
    {
        rowlen = rand() % maxrow;
        for (size_t j(0); j <= rowlen; j++)
        {
            val = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);  //all values positive
            colind = rand() % N;
            mat1.sequential_fill(colind, val);

            res[i]+=rhs.get_global(colind)*val;
        }
        mat1.end_of_row();
    }
    if (!mat1.is_filled()) LOG_ERROR("filling of right hand side and/or random Matrix failed");
    mat1.mult_vec(rhs, rhs);
    double checktol;
    checktol = std::numeric_limits<double>::epsilon() * maxrow;
    for (size_t i(fron); i <= lron; i++)
    {
        if (std::abs(rhs.get_global(i) - res[i]) > checktol)
            LOG_ERROR("spmv failed, result: ", rhs.get_global(i),
                      " ;  reference value: ", res[i], "  ; difference: ",
                       std::abs(rhs.get_global(i) - res[i]));
    }
    return 0;
}

int main()
{
    return distEllSpmv_test(100000, 30);
}
