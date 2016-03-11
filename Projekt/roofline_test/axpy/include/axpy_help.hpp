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
void axpy_check_result_(type *result, type *veconeh, type *vectwoh, size_t dim_local)
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
    for (int k = 0; k < dim_local; k++)
    {
        value += veconeh[k] * vectwoh[k];
    }
    //==========================================//
    //diff needs to be small, pretty big atm
    //==========================================//
    diff = value - result[0];
    printf(GREEN "DIFF: %f\n" RESET, diff);
    if (diff > 1.0e-2 || diff < -1.0e-2)
    {
        check = false;
    }

    if (check)
    {
        printf(GREEN "%c_Kernel outcome true\n" RESET, a);
    }
    else printf(WARNING "%c_Kernel outcome false\n", a);
}




