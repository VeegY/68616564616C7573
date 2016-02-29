#include <mpi.h>
#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "include/benchmark_help.hpp"
#include "include/timer.hpp"
using namespace std;
enum methods { unified_vs_zero, kernel_vs_kernel, kernel_vs_cpu };  //choose a method
enum version { kernel_standart, kernel_shared, kernel_advanced };   //keep your kernels in the same order as they are in the switch in your gpu_xx_call
enum memory_opt { unified, zero };                                  //choose a method of memory usage(for method k_vs_k and k_vs_cpu)
//================================================================================================/
//									GLOBAL SETTINGS!
//================================================================================================/
int method = kernel_vs_kernel;               
int version_first = kernel_standart;        
//F�r Kernel vs Kernel
int version_second = kernel_standart;      
//F�r Kernel vs Kernel und Kernel vs CPU
int memory_option = zero;

//------------------------------------------------------------------------------------------------/
//                                   APPLICATION SETTINGS
//------------------------------------------------------------------------------------------------/
#define dimlocal 16384
#define dimfvec 16384
#define maxrowlength 7
#define iteration 1000

void print_p();

template <typename type>
void cpu_ax(type *data, type *fvec, type *result, int *indices, int maxrow, int dim_local, int dim_fvec);

template<typename type>
void performance(int max_row_length, int dim_local, float time_ku, float time_ou, float time_kz, float time_oz, int runs, type schalter, int meth, int ver_first, int ver_second, int mem_option);

template<typename Scalar>
void gpu_ax_overall(Scalar *data, Scalar *fvec, Scalar *result, int *indices, int max_row_length, int dim_local, int dim_fvec, int version, int mem_option);

template<typename Scalar>
float gpu_ax_time(Scalar *data, Scalar *fvec, Scalar *result, int *indices, int max_row_length, int dim_local, int dim_fvec, int runs, int version, int mem_option);

template<typename Scalar>
void allocation(Scalar **data, Scalar **fvec, Scalar **result, int **indices, int max_row_length, int dim_local, int dim_fvec, int mem_option);

template <typename Scalar>
void cleanup(Scalar *data, Scalar *fvec, Scalar *result, int *indices, int method);

int main(int argc, char* argv[])
{
    //Array zur Zeitmessung
    //float *benchmark_times = new float[6];
    //Generiere Data/Indices Int-Array sowie fvec Int Array
    float *data_host = new float[dimlocal*maxrowlength];
    int *indices_host = new int[dimlocal*maxrowlength];
    float *fvec_host = new float[dimfvec];

    diagonal_float(data_host, indices_host, fvec_host, maxrowlength, dimlocal, dimfvec);

    Timer timer_overall,timer_cpu;

//================================================================================================/
//									THE MAGIC HAPPENS HERE
//================================================================================================/
//------------------------------------------------------------------------------------------------/
//                                Zeitmessung Overall Teil 1
//------------------------------------------------------------------------------------------------/
    if (method == unified_vs_zero) { memory_option = unified; }

    timer_overall.start();
    for(int r = 0;r<iteration;r++)
    {
        float *data_first = NULL;
        float *fvec_first = NULL;
        float *result_first = NULL;
        int *indices_first = NULL;

        allocation(&data_first, &fvec_first, &result_first, &indices_first, maxrowlength, dimlocal, dimfvec, memory_option);
        set_values(data_host, indices_host, fvec_host, data_first, indices_first, fvec_first, maxrowlength, dimlocal, dimfvec);
        gpu_ax_overall(data_first, fvec_first, result_first, indices_first, maxrowlength, dimlocal, dimfvec, version_first, memory_option);
        cleanup(data_first, fvec_first, result_first, indices_first, memory_option);
    }
    float elapsed_first_overall = timer_overall.stop() / (float)iteration;

//------------------------------------------------------------------------------------------------/
//                                Zeitmessung Overall Teil 2
//------------------------------------------------------------------------------------------------/
    if (method == unified_vs_zero) { memory_option = zero; version_second = version_first; }

    timer_overall.start();
    for (int r = 0; r<iteration; r++)
    {
        float *data_second = NULL;
        float *fvec_second = NULL;
        float *result_second = NULL;
        int *indices_second = NULL;

        if (method != kernel_vs_cpu)
        {
            allocation(&data_second, &fvec_second, &result_second, &indices_second, maxrowlength, dimlocal, dimfvec, memory_option);
            set_values(data_host, indices_host, fvec_host, data_second, indices_second, fvec_second, maxrowlength, dimlocal, dimfvec);
            gpu_ax_overall(data_second, fvec_second, result_second, indices_second, maxrowlength, dimlocal, dimfvec, version_second, memory_option);
            cleanup(data_second, fvec_second, result_second, indices_second, memory_option);
        }
        else//CPU Zeitmessung
        {
            set_values(data_host, indices_host, fvec_host, data_second, indices_second, fvec_second, maxrowlength, dimlocal, dimfvec);
            cpu_ax(data_second, fvec_second, result_second, indices_second, maxrowlength, dimlocal, dimfvec);
            cleanup(data_second, fvec_second, result_second, indices_second, 2);
        }
    }
    float elapsed_second_overall = timer_overall.stop() / (float)iteration;


//------------------------------------------------------------------------------------------------/
//                                Zeitmessung Kernel Teil 1                   
//------------------------------------------------------------------------------------------------/
    if (method == unified_vs_zero) { memory_option = unified; }
    
    float *data_first = NULL;
    float *fvec_first = NULL;
    float *result_first = NULL;
    int *indices_first = NULL;

    allocation(&data_first, &fvec_first, &result_first, &indices_first, maxrowlength, dimlocal, dimfvec, memory_option);
    set_values(data_host, indices_host, fvec_host, data_first, indices_first, fvec_first, maxrowlength, dimlocal, dimfvec);

    //=========================================//Hier muss vielleicht die Zeitmessung innerhalb der aufgerufenen Funktion stattfinden
    float elapsed_first_kernel = 
        gpu_ax_time(data_first, fvec_first, result_first, indices_first, maxrowlength, dimlocal, dimfvec, iteration, version_first, memory_option);
    //=========================================//

    cleanup(data_first, fvec_first, result_first, indices_first, memory_option);
 
 //------------------------------------------------------------------------------------------------/
 //                                Zeitmessung Kernel Teil 2                   
 //------------------------------------------------------------------------------------------------/
    if (method == unified_vs_zero) { memory_option = zero; version_second = version_first; }
    
    float *data_second = NULL;
    float *fvec_second = NULL;
    float *result_second = NULL;
    int *indices_second = NULL;
    float elapsed_second_kernel = 0.0;

    if (method != kernel_vs_cpu)
    {
        allocation(&data_second, &fvec_second, &result_second, &indices_second, maxrowlength, dimlocal, dimfvec, memory_option);
        set_values(data_host, indices_host, fvec_host, data_second, indices_second, fvec_second, maxrowlength, dimlocal, dimfvec);
        
        //=========================================//Hier muss vielleicht die Zeitmessung innerhalb der aufgerufenen Funktion stattfinden
        elapsed_second_kernel =
            gpu_ax_time(data_second, fvec_second, result_second, indices_second, maxrowlength, dimlocal, dimfvec, iteration, version_second, memory_option);
        //=========================================//
        
        cleanup(data_second, fvec_second, result_second, indices_second, memory_option);
    }
    else//CPU Zeitmessung
    {
        set_values(data_host, indices_host, fvec_host, data_second, indices_second, fvec_second, maxrowlength, dimlocal, dimfvec);
        
        //=========================================//
        timer_cpu.start();
        for (int r = 0; r < iteration; r++)
        {
            cpu_ax(data_second, fvec_second, result_second, indices_second, maxrowlength, dimlocal, dimfvec);
        }
        elapsed_second_kernel = timer_cpu.stop()*1.0e3;
        //=========================================//

        cleanup(data_second, fvec_second, result_second, indices_second, 2);
    }
    
   
//================================================================================================/
//                                         Evaluieren
//================================================================================================/

    float schalter = 0.0;
    performance(maxrowlength, dimlocal, elapsed_first_kernel, elapsed_first_overall, elapsed_second_kernel, elapsed_second_overall, iteration, schalter, method, version_first, version_second, memory_option);
  
    delete[] data_host;
    delete[] indices_host;
    delete[] fvec_host;
    return 0;
}

template <typename type>
void cpu_ax(type *data, type *fvec, type *result, int *indices, int max_row_length, int dim_local, int dim_fvec)
{
    for (int i = 0; i < dim_local; i++)
    {
        float value = 0.0;
        for (int j = 0; j < max_row_length; j++)
        {
            value += data[i + dim_local*j] * fvec[indices[i + dim_local*j]];
        }
        result[i] = value;
    }

}
template void cpu_ax<int>(int *data, int *fvec, int *result, int *indices, int max_row_length, int dim_local, int dim_fvec);
template void cpu_ax<float>(float *data, float *fvec, float *result, int *indices, int max_row_length, int dim_local, int dim_fvec);
template void cpu_ax<double>(double *data, double *fvec, double *result, int *indices, int max_row_length, int dim_local, int dim_fvec);