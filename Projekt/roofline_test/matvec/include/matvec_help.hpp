#include <stdio.h>
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <string>
#define WARNING "\e[30;41;1m"
#define GREEN "\e[32;1m"
#define RESET "\e[0m"


template <typename type>
void matvec_check_result(type *result,type *vector_host, type *data_host, int *indices_host, size_t dim, size_t max_row_length)
{
    //==========================================//
    // value is result of CPU function
    // diff is difference between GPU & CPU result
    //==========================================//
    type diff, maxdiff = 0.0;
    bool check = true;
    //==========================================//
    for (int i = 0; i < dim; i++)
    {
        type value = 0.0;
        for (int j = 0; j < max_row_length; j++)
        {
            value += data_host[i + dim*j] * vector_host[indices_host[i + dim*j]];
        }

        diff = value - result[i];
        if (diff > 1.0e-6 || diff < -1.0e-6)
        {
            check = false;
        }
        if (diff <= 0)
        {
            diff = -diff;
        }
        if (maxdiff <= diff)
        {
            maxdiff = diff;
        }
    }
    //==========================================//
    printf(GREEN "MAX.DIFF: %f\n" RESET, maxdiff);
    if (check)
    {
        printf(GREEN "Kernel outcome true\n" RESET);
    }
    else printf(WARNING "Kernel outcome false\n" RESET);
}

template<typename type>
void ellpack_fill_seven_diagonals(type *data, int *indices, int max_row_length, int dim)
{
    srand(static_cast <unsigned> (time(0)));
    int diag[7];
    diag[0] = -floor(pow(dim, (2.0 / 3.0)));
    diag[1] = -floor(pow(dim, (1.0 / 3.0)));
    diag[2] = -1;
    diag[3] = 0;
    diag[4] = 1;
    diag[5] = floor(pow(dim, (1.0 / 3.0)));
    diag[6] = floor(pow(dim, (2.0 / 3.0)));

    for (int i = 0; i < dim; i++)
    {
        int offset = 0;
        for (int j = 0; j < max_row_length; j++)
        {
            type value = 0;
            if (diag[j] >= 0 && diag[j] < dim)
            {
                value = static_cast <type> (rand()) / (static_cast <type> (RAND_MAX / 100));
                data[i + offset*dim] = value;
                indices[i + offset*dim] = diag[j];
                offset++;
            }
            diag[j] = diag[j] + 1;

        }
        for (int off = offset; off < max_row_length; off++)
        {
            data[i + offset*dim] = 0;
            indices[i + offset*dim] = 0;
        }
    }
}


