#include <mpi.h>
#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include "include/benchmark_help.hpp"
#include "include/timer.hpp"
using namespace std;
#define dim 4

template<typename Scalar>
void alloc_unified(Scalar **data, Scalar **fvec, Scalar **result, int **indices, int max_row_length, int dim_local, int dim_fvec);

template<typename Scalar>
void alloc_zero(Scalar **data, Scalar **fvec, Scalar **result, int **indices, int max_row_length, int dim_local, int dim_fvec);

template<typename Scalar>
void mult_vec_unified(Scalar* data, Scalar* fvec, Scalar* result, int* indices, int max_row_length, int dim_local, int dim_fvec);

template<typename Scalar>
void mult_vec_zero(Scalar* data, Scalar* fvec, Scalar* result, int* inices, int max_row_length, int dim_local, int dim_fvec);


int main(int argc, char* argv[])
{

//Generiere Data/Indices Int-Array sowie fvec Int Array
	int *data_host = new int[dim*dim];
	int *indices_host = new int[dim*dim];
	int *fvec_host = new int[dim];

	random_ints(data_host, indices_host, fvec_host, dim);
//Unified INT Kernel
	Timer timer_unified;

	int *data_unified = NULL;
    int *fvec_unified = NULL;
    int *result_unified = NULL;
	int *indices_unified = NULL;
	
	alloc_unified(&data_unified, &fvec_unified, &result_unified, &indices_unified, dim, dim, dim);
	set_values(data_host,indices_host,fvec_host,data_unified,indices_unified,fvec_unified, dim);
	print_stuff(data_unified, indices_unified, fvec_unified, dim);
    mult_vec_unified(data_unified, fvec_unified, result_unified, indices_unified, dim, dim, dim);

    for(int i=0;i<dim;i++)
    {
	  printf("unified float result in row %i is %f.\n", i,result_unified[i]);
	}

 
    return 0;
}
