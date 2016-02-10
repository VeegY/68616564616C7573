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
#define dimlocal 1024
#define dimfvec 1024
#define maxrowlength 7
#define iteration 1000

void print_p();

template<typename Scalar>
void performance(int max_row_length, int dim_local, float time_ku, float time_ou, float time_kz, float time_oz, int runs, Scalar schalter);

template<typename Scalar>
void alloc_unified(Scalar **data, Scalar **fvec, Scalar **result, int **indices, int max_row_length, int dim_local, int dim_fvec);

template<typename Scalar>
void alloc_zero(Scalar **data, Scalar **fvec, Scalar **result, int **indices, int max_row_length, int dim_local, int dim_fvec);

template<typename Scalar>
float mult_vec_unified_time(Scalar* data, Scalar* fvec, Scalar* result, int* indices, int max_row_length, int dim_local, int dim_fvec, int runs);

template<typename Scalar>
float mult_vec_zero_time(Scalar* data, Scalar* fvec, Scalar* result, int* inices, int max_row_length, int dim_local, int dim_fvec, int runs);

template<typename Scalar>
void mult_vec_unified(Scalar* data, Scalar* fvec, Scalar* result, int* indices, int max_row_length, int dim_local, int dim_fvec);

template<typename Scalar>
void mult_vec_zero(Scalar* data, Scalar* fvec, Scalar* result, int* inices, int max_row_length, int dim_local, int dim_fvec);

template <typename Scalar>
void cleanup(Scalar *data, Scalar *fvec, Scalar *result, int *indices, int method);

int main(int argc, char* argv[])
{
    //Generiere Data/Indices Int-Array sowie fvec Int Array
    float *data_host = new float[dimlocal*maxrowlength];
    int *indices_host = new int[dimlocal*maxrowlength];
    float *fvec_host = new float[dimfvec];

    diagonal_float(data_host, indices_host, fvec_host, maxrowlength, dimlocal, dimfvec);

    Timer timer_overall;

//================================================================================================/
//										Unified Kernel
//================================================================================================/
//------------------------------------------------------------------------------------------------/
//                                   Overall - Zeitmessung
//------------------------------------------------------------------------------------------------/
   
    timer_overall.start();
    for(int r = 0;r<iteration;r++)
    {
        float *data_unified = NULL;
        float *fvec_unified = NULL;
        float *result_unified = NULL;
        int *indices_unified = NULL;
        
        alloc_unified(&data_unified, &fvec_unified, &result_unified, &indices_unified, maxrowlength, dimlocal, dimfvec);
        set_values(data_host, indices_host, fvec_host, data_unified, indices_unified, fvec_unified, maxrowlength, dimlocal, dimfvec);
        
        mult_vec_unified(data_unified, fvec_unified, result_unified, indices_unified, maxrowlength, dimlocal, dimfvec);

        //TODO: test (0=CudaFree,1=CudeFreeHos,2=delete[])
        cleanup(data_unified, fvec_unified, result_unified, indices_unified, 0);
        //cleanup(data_unified, fvec_unified, result_unified, indices_unified, 1);
        //cleanup(data_unified, fvec_unified, result_unified, indices_unified, 2);
    }
    float elapsed_unified_overall = timer_overall.stop() / (float)iteration;

//------------------------------------------------------------------------------------------------/
//                                   Kernel - Zeitmessung
//------------------------------------------------------------------------------------------------/
    float *data_unified = NULL;
    float *fvec_unified = NULL;
    float *result_unified = NULL;
    int *indices_unified = NULL;

    alloc_unified(&data_unified, &fvec_unified, &result_unified, &indices_unified, maxrowlength, dimlocal, dimfvec);
    set_values(data_host, indices_host, fvec_host, data_unified, indices_unified, fvec_unified, maxrowlength, dimlocal, dimfvec);

    float elapsed_unified_kernel =
        mult_vec_unified_time(data_unified, fvec_unified, result_unified, indices_unified, maxrowlength, dimlocal, dimfvec, iteration);
 
    check_result(result_unified, data_host, indices_host, fvec_host, maxrowlength, dimlocal, 'u');

    //TODO: test (0=CudaFree,1=CudeFreeHos,2=delete[])
    cleanup(data_unified, fvec_unified, result_unified, indices_unified, 0);
    //cleanup(data_unified, fvec_unified, result_unified, indices_unified, 1);
    //cleanup(data_unified, fvec_unified, result_unified, indices_unified, 2);
    

//================================================================================================/
//										Zero Copy Kernel
//================================================================================================/
//------------------------------------------------------------------------------------------------/
//                                   Overall - Zeitmessung
//------------------------------------------------------------------------------------------------/

    timer_overall.start();
    for (int r = 0; r<iteration; r++)
    {
        float *data_zero = NULL;
        float *fvec_zero = NULL;
        float *result_zero = NULL;
        int *indices_zero = NULL;

        alloc_zero(&data_zero, &fvec_zero, &result_zero, &indices_zero, maxrowlength, dimlocal, dimfvec);
        set_values(data_host, indices_host, fvec_host, data_zero, indices_zero, fvec_zero, maxrowlength, dimlocal, dimfvec);

        mult_vec_zero(data_zero, fvec_zero, result_zero, indices_zero, maxrowlength, dimlocal, dimfvec);

        //TODO: test (0=CudaFree,1=CudeFreeHos,2=delete[])
        //cleanup(data_zero, fvec_zero, result_zero, indices_zero, 0);
        cleanup(data_zero, fvec_zero, result_zero, indices_zero, 1);
        //cleanup(data_zero, fvec_zero, result_zero, indices_zero, 2);
    }
    float elapsed_zero_overall = timer_overall.stop()/(float) iteration;

//------------------------------------------------------------------------------------------------/
//                                   Kernel - Zeitmessung
//------------------------------------------------------------------------------------------------/
    float *data_zero= NULL;
    float *fvec_zero = NULL;
    float *result_zero = NULL;
    int *indices_zero = NULL;

    alloc_zero(&data_zero, &fvec_zero, &result_zero, &indices_zero, maxrowlength, dimlocal, dimfvec);
    set_values(data_host, indices_host, fvec_host, data_zero, indices_zero, fvec_zero, maxrowlength, dimlocal, dimfvec);

    float elapsed_zero_kernel =
        mult_vec_zero_time(data_zero, fvec_zero, result_zero, indices_zero, maxrowlength, dimlocal, dimfvec, iteration);

    check_result(result_zero, data_host, indices_host, fvec_host, maxrowlength, dimlocal, 'z');

    //TODO: test (0=CudaFree,1=CudeFreeHos,2=delete[])
    //cleanup(data_zero, fvec_zero, result_zero, indices_zero, 0);
    cleanup(data_zero, fvec_zero, result_zero, indices_zero, 1);
    //cleanup(data_zero, fvec_zero, result_zero, indices_zero, 2);

//================================================================================================/
//                                         Evaluieren
//================================================================================================/
    
    //print_p();
    float schalter = 0.0;
    performance(maxrowlength, dimlocal, elapsed_unified_kernel, elapsed_unified_overall, elapsed_zero_kernel, elapsed_zero_overall, iteration, schalter);

    delete[] data_host;
    delete[] indices_host;
    delete[] fvec_host;
    return 0;
}
