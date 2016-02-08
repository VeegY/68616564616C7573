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
#define dimlocal 8
#define dimfvec 8
#define maxrowlength 7
#define iteration 1

void print_p();

template<typename Scalar>
void performance(int max_row_length, int dim_local, float time_ku, float time_ou, Scalar schalter);

template<typename Scalar>
void alloc_unified(Scalar **data, Scalar **fvec, Scalar **result, int **indices, int max_row_length, int dim_local, int dim_fvec);

template<typename Scalar>
void alloc_zero(Scalar **data, Scalar **fvec, Scalar **result, int **indices, int max_row_length, int dim_local, int dim_fvec);

template<typename Scalar>
void mult_vec_unified_time(Scalar* data, Scalar* fvec, Scalar* result, int* indices, int max_row_length, int dim_local, int dim_fvec, int runs);

template<typename Scalar>
void mult_vec_zero_time(Scalar* data, Scalar* fvec, Scalar* result, int* inices, int max_row_length, int dim_local, int dim_fvec, int runs);

template<typename Scalar>
void mult_vec_unified(Scalar* data, Scalar* fvec, Scalar* result, int* indices, int max_row_length, int dim_local, int dim_fvec);

template<typename Scalar>
void mult_vec_zero(Scalar* data, Scalar* fvec, Scalar* result, int* inices, int max_row_length, int dim_local, int dim_fvec);

template <typename Scalar>
void cleanup(Scalar *data, Scalar *fvec, Scalar *result, int *indices, int method);

int main(int argc, char* argv[])
{
    cout << "DIM = " << dimlocal << " - RUNS = " << iteration << "\n";

//Generiere Data/Indices Int-Array sowie fvec Int Array
    float *data_host = new float[dimlocal*maxrowlength];
    int *indices_host = new int[dimlocal*maxrowlength];
    float *fvec_host = new float[dimfvec];
    float *result_host = new float[dimlocal];

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
    float elapsed_unified_overall = timer_overall.stop();

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

    check_result(result_unified, data_host, indices_host, fvec_host, maxrowlength, dimlocal, "U");

    //TODO: test (0=CudaFree,1=CudeFreeHos,2=delete[])
    cleanup(data_unified, fvec_unified, result_unified, indices_unified, 0);
    //cleanup(data_unified, fvec_unified, result_unified, indices_unified, 1);
    //cleanup(data_unified, fvec_unified, result_unified, indices_unified, 2);
    

//================================================================================================/
//										Zero Copy Kernel
//================================================================================================/
/*for(int r = 0;r<runs;r++)
{
    //Initialisiere Timer und Pointer
    timer_overall.start();

    int *data_zero = NULL;
    int *fvec_zero = NULL;
    int *result_zero = NULL;
    int *indices_zero = NULL;


    //Setze Pointer als ZERO_COPY Pointer fest und fuelle mit Daten
    alloc_zero(&data_zero, &fvec_zero, &result_zero, &indices_zero, dim, dim, dim);
    set_values(data_host, indices_host, fvec_host, data_zero, indices_zero, fvec_zero, dim);


    //Starte Kernel Stoppuhr und Kernel
    timer_kernel.start();
    mult_vec_zero(data_zero, fvec_zero, result_zero, indices_zero, dim, dim, dim);


    //Messe Zeiten
    zero_kernel_time[r] = timer_kernel.stop();
    zero_overall_time[r] = timer_overall.stop();


    //Alte Zeitmessung
    //float elapsed_zero_kernel = timer_zero_kernel.stop();
    //float elapsed_zero_overall = timer_zero_overall.stop();
    //cout << "KERNEL TIME: " << elapsed_zero_kernel * 1000 << "\n";
    //cout << "OVERALL TIMER: " << elapsed_zero_overall * 1000 << "\n\n";


    //Ein Ergebnis Check insgesamt
    if(r==0)
    {
      if(check_result(result_zero, data_host, indices_host, fvec_host, dim))
      {
          cout << "*Z_CORRECT*\n";
      }
          else{cout << "*Z_FALSE*\n";}
    }


    //Aufraeumen
    cleanup(data_zero, fvec_zero, result_zero, indices_zero);

}
//================================================================================================/*/

    //print_time(unified_kernel_time, unified_overall_time,zero_kernel_time,zero_overall_time,runs);


    delete[] data_host;
    delete[] indices_host;
    delete[] fvec_host;
    return 0;
}
