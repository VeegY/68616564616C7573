#include <stdio.h>
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <string>
#define GREEN "\e[32;1m"
#define RESET "\e[0m"
#define BLUE "\e[34;1m"
#define CYAN "\e[36;1m"
#define GREY "\e[30;1m"
#define MAGENTA "\e[35;1m"

template<typename type>
void performance(int dim, float overall, float kernel, int runs, type schalter)
{
    int bytes = 2 * sizeof(type) * dim;
    int flop = dim;
    double ai = (1 / 16);

    printf(GREY    "===============================================\n");
    printf(MAGENTA "                PERFORMANCE\n");
    printf(        "                COPY KERNEL\n");
    printf(GREY    "===============================================\n");
    printf(        "-----------------------------------------------\n");
    printf(GREY    "-----------------------------------------------\n");
    printf(CYAN    "Kernel Runtime:\t\t\t%f(ms)\n", kernel);
    printf(        "Overall Runtime:\t\t%f(ms)\n", overall);
    printf(        "Bandwith(th. Peak):\t\t%.2f(14.9)(GB/s)\n", bytes / (time_ku*1.0e6));
    printf(        "Flops(th. Peak):\t\t%.6f(326)(GFLOPS/s)\n", flop / (time_ku*1.0e6));
    printf(        "StreamBW * AI:\t\t\t%f\n", 50 * (ai));
    printf(GREY    "-----------------------------------------------\n");
    printf(        "-----------------------------------------------\n");

}
template void performance<int>(int dim, float overall, float kernel, int runs, int schalter);
template void performance<float>(int dim, float overall, float kernel, int runs, float schalter);
template void performance<double>(int dim, float overall, float kernel, int runs, double schalter);


//=============================================================================
///////////////////////////////////////////////////////////////////////////////
///                             HELPER FUNCTIONS                            ///
///////////////////////////////////////////////////////////////////////////////                       
//=============================================================================
template <typename type>
void set_data(type *data, size_t size)
{
    srand(static_cast <unsigned> (time(0)));
    type value;
    for (int k = 0; k < size; k++)
    {
        value = static_cast <type> (rand()) / (static_cast <type> (RAND_MAX / 10));
        data[k] = value;
    }

}

template <typename type>
void copy_data(type *h_data, type *d_data, size_t size)
{
    for (int k = 0; k < size; k++)
    {
        d_data[k] = h_data[k];
    }
}

template <typename type>
void print_data(type *data, std::string str, size_t size)
{
    using std::string;
    for (int k = 0; k < size; k++)
    {
        printf("%s\t%i\t%d\n", str.c.str(), k, data[k]);
    }
}