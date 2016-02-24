#include <mpi.h>
#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "include/profile_help.hpp"
#include "include/timer.hpp"
using namespace std;
#define dimlocal 16384
#define dimfvec 16384
#define maxrowlength 7
#define iteration 1

void print_p();

template<typename Scalar>
void performance(int max_row_length, int dim_local, float time_ku, float time_ou, float time_kz, float time_oz, int runs, Scalar schalter);


//======================ALLOCATION UNIFIED MEMORY==================//
template<typename Scalar>
void alloc_unified(Scalar **data, Scalar **fvec, Scalar **result, int **indices, int max_row_length, int dim_local, int dim_fvec);

template<typename Scalar>
float profile_alloc_unified(Scalar **data, Scalar **fvec, Scalar **result, int **indices, int max_row_length, int dim_local, int dim_fvec);


//======================ALLOCATION ZERO COPY=======================//
template<typename Scalar>
void alloc_zero(Scalar **data, Scalar **fvec, Scalar **result, int **indices, int max_row_length, int dim_local, int dim_fvec);

template<typename Scalar>
float profile_alloc_zero(Scalar **data, Scalar **fvec, Scalar **result, int ** indices, int max_row_length, int dim_local, int dim_fvec);


//======================UNIFIED MEMORY KERNEL CALLS================//
template<typename Scalar>
void mult_vec_unified(Scalar* data, Scalar* fvec, Scalar* result, int* indices, int max_row_length, int dim_local, int dim_fvec);

template<typename Scalar>
float mult_vec_unified_time(Scalar* data, Scalar* fvec, Scalar* result, int* indices, int max_row_length, int dim_local, int dim_fvec, int runs);


//======================ZERO COPY KERNEL CALLS=====================//
template<typename Scalar>
void mult_vec_zero(Scalar* data, Scalar* fvec, Scalar* result, int* inices, int max_row_length, int dim_local, int dim_fvec);

template<typename Scalar>
float mult_vec_zero_time(Scalar* data, Scalar* fvec, Scalar* result, int* inices, int max_row_length, int dim_local, int dim_fvec, int runs);

template<typename Scalar>
void profile_mult_vec_zero(Scalar *data, Scalar *fvec, Scalar *result, int *indices, int max_row_length, int dim_local, int dim_fvec, int runs, float *profile);


//======================CLEANUP====================================//
template <typename Scalar>
void cleanup(Scalar *data, Scalar *fvec, Scalar *result, int *indices, int method);

template <typename Scalar>
float profile_cleanup(Scalar *data, Scalar *fvec, Scalar *result, int *indices, int method)



int main(int argc, char* argv[])
{
    //Generiere Data/Indices Int-Array sowie fvec Int Array
    float *data_host = new float[dimlocal*maxrowlength];
    int *indices_host = new int[dimlocal*maxrowlength];
    float *fvec_host = new float[dimfvec];

    diagonal_float(data_host, indices_host, fvec_host, maxrowlength, dimlocal, dimfvec);

    /*Timer timer_overall;

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

        //TODO: test (0=CudaFree,1=CudeFreeHost,2=delete[])
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
    float elapsed_unified_kernel = 1.0, elapsed_unified_overall = 1.0;
    print_p();
    float schalter = 0.0;
    performance(maxrowlength, dimlocal, elapsed_unified_kernel, elapsed_unified_overall, elapsed_zero_kernel, elapsed_zero_overall, iteration, schalter);
    */

    //================================================================================================/
    //                                  Profile
    //================================================================================================/

    Timer set_val_unified, set_val_zero;
    float *profile = new profile[10];
    for (int k = 0; k < 10; k++)
    {
        profile[k] = 0.0;
    }

    //------------------------------------------------------------------------------------------------/
    //                                  Unified Memory
    //------------------------------------------------------------------------------------------------/
    float *data_unified = NULL;
    float *fvec_unified = NULL;
    float *result_unified = NULL;
    int *indices_unified = NULL;

    //===============================//PROFILE TIME [0] ~ allocation
    profile[0] += profile_alloc_unified(&data_unified, &fvec_unified, &result_unified, &indices_unified, maxrowlength, dimlocal, dimfvec);
    //===============================//

    //===============================//PROFILE TIME [1] ~ set_values
    set_val_unified.start();
    set_values(data_host, indices_host, fvec_host, data_unified, indices_unified, fvec_unified, maxrowlength, dimlocal, dimfvec);
    profile[1] += set_val_unified.stop();
    //===============================//

    //===============================//PROFILE TIME [2] ~ kernel
    profile[2] += mult_vec_unified(data_unified, fvec_unified, result_unified, indices_unified, maxrowlength, dimlocal, dimfvec, iteration);
    //===============================//


    check_result(result_unified, data_host, indices_host, fvec_host, maxrowlength, dimlocal, 'u');

    //===============================//PROFILE TIME [3] ~ cleanup
    //TODO: test (0=CudaFree,1=CudeFreeHos,2=delete[])
    profile[3] += profile_cleanup(data_unified, fvec_unified, result_unified, indices_unified, 0);
    //cleanup(data_unified, fvec_unified, result_unified, indices_unified, 1);
    //cleanup(data_unified, fvec_unified, result_unified, indices_unified, 2);


//------------------------------------------------------------------------------------------------/
//                                  Zero Copy
//------------------------------------------------------------------------------------------------/

    float *data_zero = NULL;
    float *fvec_zero = NULL;
    float *result_zero = NULL;
    int *indices_zero = NULL;

    //===============================//PROFILE TIME [4] ~ allocate
    profile[4] += profile_alloc_zero(&data_zero, &fvec_zero, &result_zero, &indices_zero, maxrowlength, dimlocal, dimfvec);
    //===============================//

    //===============================//PROFILE TIME [5] ~ set_values
    set_val_zero.start();
    set_values(data_host, indices_host, fvec_host, data_zero, indices_zero, fvec_zero, maxrowlength, dimlocal, dimfvec);
    profile[5] += set_val_zero.stop()*1.0e3;
    //===============================//

    //===============================//PROFILE TIME [6] - PROFILE TIME [8] ~ gethostpointer ~ kernel ~ cleanup
    profile_mult_vec_zero(data_zero, fvec_zero, result_zero, indices_zero, maxrowlength, dimlocal, dimfvec, iteration, profile);
    //===============================//

    check_result(result_zero, data_host, indices_host, fvec_host, maxrowlength, dimlocal, 'z');

    //===============================//PROFILE TIME [9] ~ cleanup
    //TODO: test (0=CudaFree,1=CudeFreeHos,2=delete[])
    //cleanup(data_zero, fvec_zero, result_zero, indices_zero, 0);
    profile[9] += profile_cleanup(data_zero, fvec_zero, result_zero, indices_zero, 1);
    //cleanup(data_zero, fvec_zero, result_zero, indices_zero, 2);
    //===============================//
    
    for (int i = 0; i < 10; i++)
    {
        printf("profile %i: %fms\n", i, profile[i]/iteration);
    }
    
    
    
    delete[] data_host;
    delete[] indices_host;
    delete[] fvec_host;
    return 0;
}
