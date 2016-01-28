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
#define dim 64
#define runs 1

void print_p();

template<typename Scalar>
void alloc_unified(Scalar **data, Scalar **fvec, Scalar **result, int **indices, int max_row_length, int dim_local, int dim_fvec);

template<typename Scalar>
void alloc_zero(Scalar **data, Scalar **fvec, Scalar **result, int **indices, int max_row_length, int dim_local, int dim_fvec);

template<typename Scalar>
void mult_vec_unified(Scalar* data, Scalar* fvec, Scalar* result, int* indices, int max_row_length, int dim_local, int dim_fvec);

template<typename Scalar>
void mult_vec_zero(Scalar* data, Scalar* fvec, Scalar* result, int* inices, int max_row_length, int dim_local, int dim_fvec);

template <typename Scalar>
void cleanup(Scalar *data, Scalar *fvec, Scalar *result, int *indices);

int main(int argc, char* argv[])
{

    cout << "DIM = " << dim << " - RUNS = " << runs << "\n";
    float *unified_kernel_time = new float[runs];
    float *unified_overall_time = new float[runs];
    float *zero_kernel_time = new float[runs];
    float *zero_overall_time = new float[runs];

//Generiere Data/Indices Int-Array sowie fvec Int Array
    int *data_host = new int[dim*dim];
    int *indices_host = new int[dim*dim];
    int *fvec_host = new int[dim];

    random_ints(data_host, indices_host, fvec_host, dim);


    Timer timer_kernel;
    Timer timer_overall;

//================================================================================================/
//										Unified Kernel
//================================================================================================/
for(int r = 0;r<runs;r++)
{
    //Initialisiere Timer und Pointer
    timer_overall.start();

    int *data_unified = NULL;
    int *fvec_unified = NULL;
    int *result_unified = NULL;
    int *indices_unified = NULL;


    //Setze Pointer als UNIFIED Pointer mit bestimmter Groesze und setze Daten
    alloc_unified(&data_unified, &fvec_unified, &result_unified, &indices_unified, dim, dim, dim);
    set_values(data_host,indices_host,fvec_host,data_unified,indices_unified,fvec_unified, dim);


    //Kernel Zeit und Kernel starten
    timer_kernel.start();
    mult_vec_unified(data_unified, fvec_unified, result_unified, indices_unified, dim, dim, dim);


    //Messe Zeiten
    float elapsed_k = timer_kernel.stop();
    float elapsed_o = timer_overall.stop();
    unified_kernel_time[r] = elapsed_k;
    unified_overall_time[r] = elapsed_o;

    /*//TIME
    float elapsed_unified_kernel = timer_unified_kernel.stop();
    float elapsed_unfified_overall = timer_unified_overall.stop();
    cout << "KERNEL TIME: " << elapsed_unified_kernel * 1000 << "\n";
    cout << "OVERALL TIMER: " << elapsed_unfified_overall * 1000 << "\n\n";*/


    //Ein Ergebnis Check pro Schleife
    // print_stuff(data_unified, indices_unified, fvec_unified,result_unified, dim);
    if(r==0)
    {
      if(check_result(result_unified, data_host, indices_host, fvec_host, dim))
      {
          cout << "*U_CORRECT*\n";
      }
      else{cout << "*U_FALSE*\n";}
    }


    //Aufraumen und Speicher freigebene
    cleanup(data_unified, fvec_unified, result_unified, indices_unified);

}
//================================================================================================/
//										Zero Copy Kernel
//================================================================================================/
for(int r = 0;r<runs;r++)
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
//================================================================================================/

    print_time(unified_kernel_time, unified_overall_time,zero_kernel_time,zero_overall_time,runs);


    delete[] data_host;
    delete[] indices_host;
    delete[] fvec_host;
    return 0;
}
