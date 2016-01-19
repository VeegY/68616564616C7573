#include <mpi.h>
#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <iomanip>
using namespace std;

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

     MPI_Init(&argc, &argv);
     //MPI_Status status;

     int rank;
     MPI_Comm_rank(MPI_COMM_WORLD, &rank);

     if(rank == 0)
     {
       int max_row_length = 2;
       int dim_local = 2;
       int dim_fvec = 2;

       float *data = NULL;
       float *fvec = NULL;
       float *result = NULL;
       int *indices = NULL;

       alloc_unified(&data,&fvec,&result,&indices,max_row_length,dim_local,dim_fvec);

       data[0]=1;
       data[1]=2;
       data[2]=3;
       data[3]=4;

       fvec[0]=1;
       fvec[1]=2;

       indices[0]=0;
       indices[1]=1;
       indices[2]=0;
       indices[3]=1;

       mult_vec_unified(data, fvec, result, indices, max_row_length, dim_local, dim_fvec);

       for(int i=0;i<dim_fvec;i++)
       {
         printf("Rank %i unified float result in row %i is %f.\n",rank , i,result[i]);
       }

     }


     if(rank == 1)
     {
       int max_row_length = 2;
       int dim_local = 2;
       int dim_fvec = 2;

       int *data = NULL;
       int *fvec = NULL;
       int *result = NULL;
       int *indices = NULL;

       alloc_zero(&data,&fvec,&result,&indices,max_row_length,dim_local,dim_fvec);

       data[0]=1;
       data[1]=2;
       data[2]=3;
       data[3]=4;

       fvec[0]=1;
       fvec[1]=2;

       indices[0]=0;
       indices[1]=1;
       indices[2]=0;
       indices[3]=1;

       mult_vec_zero(data, fvec, result, indices, max_row_length, dim_local, dim_fvec);

       for(int i=0;i<dim_fvec;i++)
       {
         printf("Rank %i zero copy int result in row %i is %i.\n", rank, i, result[i]);
       }

     }
     MPI_Finalize();
     return 0;
}
