#include <mpi.h>
#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include "include/benchmark_help.hpp"
#include "include/timer.hpp"
using namespace std;
#define dim 3000

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

//Generiere Data/Indices Int-Array sowie fvec Int Array
    int *data_host = new int[dim*dim];
    int *indices_host = new int[dim*dim];
    int *fvec_host = new int[dim];

    random_ints(data_host, indices_host, fvec_host, dim);
    
//================================================================================================/
//										Unified Kernel
//================================================================================================/

	Timer timer_unified_kernel;
        Timer timer_unified_overall;
	cout << "UNIFIED KERNEL STARTED" << "\n--------------------------------------------------------------\n";
	timer_unified_overall.start();

    int *data_unified = NULL;
    int *fvec_unified = NULL;
    int *result_unified = NULL;
    int *indices_unified = NULL;
    
    alloc_unified(&data_unified, &fvec_unified, &result_unified, &indices_unified, dim, dim, dim);
    set_values(data_host,indices_host,fvec_host,data_unified,indices_unified,fvec_unified, dim);

    timer_unified_kernel.start();

    mult_vec_unified(data_unified, fvec_unified, result_unified, indices_unified, dim, dim, dim);

	float elapsed_unified_kernel = timer_unified_kernel.stop();
	float elapsed_unfified_overall = timer_unified_overall.stop();
    cout << "KERNEL TIME: " << elapsed_unified_kernel * 1000 << "\n";
	cout << "OVERALL TIMER: " << elapsed_unfified_overall * 1000 << "\n\n";

	//cleanup(data_unified, fvec_unified, result_unified, indices_unified);


//================================================================================================/
//										Zero Copy Kernel
//================================================================================================/

	Timer timer_zero_kernel;
        Timer timer_zero_overall;
	cout << "ZERO_COPY KERNEL STARTED" << "\n--------------------------------------------------------------\n";
	timer_zero_overall.start();

	int *data_zero = NULL;
	int *fvec_zero = NULL;
	int *result_zero = NULL;
	int *indices_zero = NULL;

	alloc_zero(&data_zero, &fvec_zero, &result_zero, &indices_zero, dim, dim, dim);
	set_values(data_host, indices_host, fvec_host, data_zero, indices_zero, fvec_zero, dim);

	timer_zero_kernel.start();

	mult_vec_zero(data_zero, fvec_zero, result_zero, indices_zero, dim, dim, dim);

	float elapsed_zero_kernel = timer_zero_kernel.stop();
	float elapsed_zero_overall = timer_zero_overall.stop();
	cout << "KERNEL TIME: " << elapsed_zero_kernel * 1000 << "\n";
	cout << "OVERALL TIMER: " << elapsed_zero_overall * 1000 << "\n\n";

	cleanup(data_zero, fvec_zero, result_zero, indices_zero);

//================================================================================================/

    delete[] data_host;
    delete[] indices_host;
    delete[] fvec_host;
    return 0;
}
