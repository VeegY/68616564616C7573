#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <string.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cmath>

//KERNEL
template<typename type>
__global__ void  gpu_ax(type* data, type* fvec, type* result, int* indices, int max_row_length, int dim_local)
{

    int idx = threadIdx.x+blockDim.x*blockIdx.x;
    if(idx<dim_local)
    {
      int col;
      type svalue = 0, value;
      for(int i = 0;i < max_row_length; i++)
      {
        value = data[i*dim_local+idx];
        col = indices[i*dim_local+idx];
        svalue += value*fvec[col];
      }
      result[idx]=svalue;
    }
}

//PROPERTIES OF TEGRA K1
void print_p()
{

    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop,0);

    printf("Max Threads Per Block: %i\n", prop.maxThreadsPerBlock);
    printf("Max Grid Size: %ix%ix%i\n", prop.maxGridSize[0],prop.maxGridSize[1],prop.maxGridSize[2]);

}

//ALLOCATE MEMORY FUNCTION FOR UNIFIED MEMORY
template<typename Scalar>
void alloc_unified(Scalar **data, Scalar **fvec, Scalar **result, int **indices, int max_row_length, int dim_local,int dim_fvec)
{
    cudaMallocManaged((void **)data, sizeof(Scalar)*dim_local*max_row_length);
    cudaMallocManaged((void **)fvec, sizeof(Scalar)*dim_fvec);
    cudaMallocManaged((void **)result, sizeof(Scalar)*dim_local);
    cudaMallocManaged((void **)indices, sizeof(int)*dim_local*max_row_length);
}
template void alloc_unified<int>(int **data, int **fvec, int **result, int **indices, int max_row_length, int dim_local, int dim_fvec);
template void alloc_unified<float>(float **data, float **fvec, float **result, int **indices, int max_row_length, int dim_local, int dim_fvec);
template void alloc_unified<double>(double **data, double **fvec, double **result, int **indices, int max_row_length, int dim_local, int dim_fvec);


//ALLOCATE MEMORY FUNCTION FOR ZERO COPY 
template<typename Scalar>
void alloc_zero(Scalar **data, Scalar **fvec, Scalar **result, int ** indices, int max_row_length, int dim_local, int dim_fvec)
{
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop,0);

    if(prop.canMapHostMemory)
    {
      cudaSetDeviceFlags(cudaDeviceMapHost);

      cudaHostAlloc((void **)data, sizeof(Scalar)*max_row_length*dim_local, cudaHostAllocMapped);
      cudaHostAlloc((void **)fvec, sizeof(Scalar)*dim_fvec, cudaHostAllocMapped);
      cudaHostAlloc((void **)result, sizeof(Scalar)*dim_local, cudaHostAllocMapped);
      cudaHostAlloc((void **)indices, sizeof(int)*max_row_length*dim_local, cudaHostAllocMapped);
    }
}
template void alloc_zero<int>(int **data, int **fvec, int **result, int **indices, int max_row_length, int dim_local, int dim_fvec);
template void alloc_zero<float>(float **data, float **fvec, float **result, int **indices, int max_row_length, int dim_local, int dim_fvec);
template void alloc_zero<double>(double **data, double **fvec, double **result, int **indices, int max_row_length, int dim_local, int dim_fvec);


//KERNEL CALL WITH UNIFIED MEMORY (NEED TO CALL ALLOC_UNIFIED BEFORE)
template<typename Scalar>
void mult_vec_unified(Scalar *data, Scalar *fvec, Scalar *result, int *indices, int max_row_length, int dim_local, int dim_fvec)
{
    int num_blocks = ceil((double)dim_local/1024);
    int num_threads = ceil(((double)dim_local/num_blocks)/32)*32;

    gpu_ax<<<num_blocks,num_threads>>>(data,fvec,result,indices,max_row_length, dim_local);
    cudaDeviceSynchronize();

}
template void mult_vec_unified<int>(int* data, int* fvec, int* result, int* indices, int max_row_length, int dim_local,int dim_fvec);
template void mult_vec_unified<float>(float* data, float* fvec, float* result, int* indices, int max_row_length, int dim_local, int dim_fvec);
template void mult_vec_unified<double>(double* data, double* fvec, double* restult, int* indices, int max_row_length, int dim_local, int dim_fvec);


//KERNE CALL WITH ZERO COPY (NEED TO CALL ALLOC_ZERO BEFORE)
template<typename Scalar>
void mult_vec_zero(Scalar *data, Scalar *fvec, Scalar *result, int *indices, int max_row_length, int dim_local, int dim_fvec)
{
    Scalar *d_data, *d_fvec, *d_result;
    int *d_indices;

    cudaHostGetDevicePointer((void **)&d_data,(void *)data, 0);
    cudaHostGetDevicePointer((void **)&d_fvec, (void *)fvec, 0);
    cudaHostGetDevicePointer((void **)&d_result, (void *)result, 0);
    cudaHostGetDevicePointer((void **)&d_indices, (void *)indices, 0);

    int num_blocks = ceil((double)dim_local/1024);
    int num_threads = ceil(((double)dim_local/num_blocks)/32)*32;

    gpu_ax<<<num_blocks,num_threads>>>(d_data, d_fvec, d_result, d_indices, max_row_length, dim_local);
    cudaDeviceSynchronize();
}
template void mult_vec_zero<int>(int* data, int* fvec, int* result, int* indices, int max_row_length, int dim_local, int  dim_fvec);
template void mult_vec_zero<float>(float* data, float* fvec, float* result, int* indices, int max_row_length, int dim_local, int dim_fvec);
template void mult_vec_zero<double>(double* data, double* fvec, double* restult, int* indices, int max_row_length, int dim_local, int dim_fvec);

template <typename Scalar>
void cleanup(Scalar *data, Scalar *fvec, Scalar *result, int *indices)
{
	cudaFreeHost(data);
	cudaFreeHost(fvec);
	cudaFreeHost(result);
	cudaFreeHost(indices);
}
template void cleanup<int>(int *data, int *fvec, int *result, int *indices);
template void cleanup<float>(float *data, float *fvec, float *result, int *indices);
template void cleanup<double>(double *data, double *fvec, double *result, int *indices);
