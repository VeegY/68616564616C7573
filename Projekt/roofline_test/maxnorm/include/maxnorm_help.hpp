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
void maxnorm_check_result_(type *result, type *vector_host, size_t dim_local)
{
    //==========================================//
    // value is result of CPU function
    // diff is difference between GPU & CPU result
    //==========================================//
    type diff, value = 0.0, current = 0.0;
    bool check = true;
    //==========================================//
    for (int i = 0; i < dim_local; i++)
    {
        current = vector[i];
        if (current <= 0)
        {
            current = -current;
        }
        if (current >= value)
        {
            value = current;
        }
    }
    //==========================================//
    diff = value) - result[0];
    printf(GREEN "DIFF: %f\n" RESET, diff);
    if (diff > 1.0e-6 || diff < -1.0e-6)
    {
        check = false;
    }

    if (check)
    {
        printf(GREEN "Kernel outcome true\n" RESET);
    }
    else printf(WARNING "Kernel outcome false\n" RESET);
}




