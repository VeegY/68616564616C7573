#include <stdio.h>
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <ctime>
#include <cstdlib>
#define GREEN "\e[32;1m"
#define RESET "\e[0m"
#define BLUE "\e[34;1m"
#define CYAN "\e[36;1m"
#define GREY "\e[30;1m"
#define MAGENTA "\e[35;1m"

template<typename type>
void performance(int max_row_length, int dim_local, float time_ku, float time_ou, float time_kz, float time_oz, int runs, type schalter, int meth, int ver_first, int ver_second, int mem_option)
{
    //===#ELEMENTS IN THE MATRIX===================================//
    unsigned long long int elements = 7 * dim_local - 2 - 2 * (floor(pow(dim_local, (1.0 / 3.0)))) - 2 * (floor(pow(dim_local, (2.0 / 3.0))));

    //==='DISK STORAGE~============================================//
    unsigned long long int storage = sizeof(type)*(2 * dim_local + dim_local*max_row_length) + sizeof(int)*dim_local*max_row_length;

    //===#FLOP=====================================================//
    unsigned long long int flop = 2 * elements;

    //==#BYTES=====================================================//           
    int bytes = elements*(sizeof(type) + sizeof(int)) + 2 * (sizeof(type)*dim_local);// Elements(Data+Indices) + Fvec Read und Result Write
    printf(GREY "===============================================\n");
    printf(MAGENTA "                PERFORMANCE\n");
    printf("           %s\n", method.c_str());
    printf("        DIM = %i ~~ %i Iterations\n", dim_local, runs);
    printf("            %.2fGB/2GB DRAM used\n", storage / 1.0e9);
    printf(GREY "===============================================\n");
    printf("-----------------------------------------------\n");
    printf(CYAN "                    %s\n", first.c_str());
    printf(GREY "-----------------------------------------------\n");
    printf(CYAN "Kernel Runtime:\t\t\t%f(ms)\n", time_ku);
    printf("Overall Runtime:\t\t%f(ms)\n", time_ou*1.0e3);
    printf("Bandwith(th. Peak):\t\t%.2f(14.9)(GB/s)\n", bytes / (time_ku*1.0e6));
    printf("Flops(th. Peak):\t\t%.6f(326)(GFLOPS/s)\n", flop / (time_ku*1.0e6));
    printf(GREY "-----------------------------------------------\n");
    printf("-----------------------------------------------\n");
    printf(BLUE "                     %s\n", second.c_str());
    printf(GREY "-----------------------------------------------\n");
    printf(BLUE "Kernel Runtime:\t\t\t%f(ms)\n", time_kz);
    printf("Overall Runtime:\t\t%f(ms)\n", time_oz*1.0e3);
    printf("Bandwith(th. Peak):\t\t%.2f(14.9)(GB/s)\n", bytes / (time_kz*1.0e6));
    printf("Flops(th. Peak):\t\t%.6f(326)(GFLOPS/s)\n", flop / (time_kz*1.0e6));
    printf(GREY "-----------------------------------------------\n" RESET);

}
template void performance<int>(int max_row_length, int dim_local, float time_ku, float time_ou, float time_kz, float time_oz, int runs, int schalter, int meth, int ver_first, int ver_second, int mem_option);
template void performance<float>(int max_row_length, int dim_local, float time_ku, float time_ou, float time_kz, float time_oz, int runs, float schalter, int meth, int ver_first, int ver_second, int mem_option);
template void performance<double>(int max_row_length, int dim_local, float time_ku, float time_ou, float time_kz, float time_oz, int runs, double schalter, int meth, int ver_first, int ver_second, int mem_option);


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
void print_data(type *data, string str, size_t size)
{
    using std::string;
    for (int k = 0; k < size; k++)
    {
        printf("%s\t%i\t%d\n", str.c.str(), k, data[k]);
    }
}