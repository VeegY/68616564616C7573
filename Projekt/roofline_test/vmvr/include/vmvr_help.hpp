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
void dotproduct_check_result(type *placehold, type *veconeh, type *vectwoh, size_t dim_local)
{
    int numblock = ceil((double)dim_local / 1024);
    //==========================================//
    // value is result of CPU function
    // diff is difference between GPU & CPU result
    //==========================================//
    type diff, maxdiff = 0.0;
    bool check = true;
    //==========================================//
    for (int i = 0; i < numblock; i++)
    {
        type value = 0.0;
        for (int j = 0; j < 1024; j++)
        {
            value += veconeh[i * 1024 + j] * vectwoh[i * 1024 + j];
        }

        diff = value - placehold[i];
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
    else printf(WARNING "Kernel outcome false\n");
}

template <typename type>
void l2norm_check_result(type *placehold, type *vector, size_t dim_local)
{
    int numblock = ceil((double)dim_local / 1024);
    //==========================================//
    // value is result of CPU function
    // diff is difference between GPU & CPU result
    //==========================================//
    type diff, maxdiff = 0.0;
    bool check = true;
    //==========================================//
    for (int i = 0; i < numblock; i++)
    {
        type value = 0.0;
        for (int j = 0; j < 1024; j++)
        {
            value += vector[i * 1024 + j] * vector[i * 1024 + j];
        }

        diff = value - placehold[i];
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
    else printf(WARNING "Kernel outcome false\n");
}

template <typename type>
void reduce_check_result(type *result, type *placehold, size_t dim_local)
{
    //==========================================//
    // value is result of CPU function
    // diff is difference between GPU & CPU result
    //==========================================//
    type diff, value = 0.0;
    bool check = true;
    //==========================================//
    //calculate CPU result
    //==========================================//
    for (int i = 0; i < dim_local; i++)
    {
        //printf("placehold %i :%f\n", i, placehold[i]);
        value += placehold[i];
    }
    //==========================================//
    //diff needs to be small, pretty big atm
    //==========================================//
    value = sqrt(value);
    diff = value - result[0];
    printf(GREEN "DIFF: %f\n" RESET, diff);
    if (diff > 1.0e-6 || diff < -1.0e-6)
    {
        check = false;
    }

    if (check)
    {
        printf(GREEN "Kernel outcome true\n" RESET);
    }
    else printf("Kernel outcome false\n");
}
