#include <seqan/parallel.h>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <string.h>
#include <cuda.h>
#include <cuda_runtime.h>
//#include "timer.hpp"
#define N 5


template<typename Scalar>
__global__ void  gpu_ax(Scalar *data, Scalar* fvec, Scalar* result, int *indices)
{
    int idx = threadIdx.x+blockIdx.x*blockDim.x;
    int row = blockIdx.x;

    Scalar value = 0;
    if(!data[idx] == 0)
    {
      value = data[idx]*fvec[indices[idx]];
      atomicAdd(&result[row],value);
    }

}


__global__ void  gpu_ax(double *data, double* fvec, double* result, int *indices)
{
    int idx = threadIdx.x+blockIdx.x*blockDim.x;
    int row = blockIdx.x;

    double value = 0;
    if(!data[idx] == 0)
    {
      value = data[idx]*fvec[indices[idx]];
      atomicCAS(&result[row],result[row],result[row]+value);
    }
}

template<typename Scalar>
void mult_vec_unified(Scalar *data, Scalar *fvec, Scalar *result, int *indices, int max_row_length, int dim_local)
{
//    Timer timer;
//    timer.start();

    cudaMallocManaged(&data, sizeof(Scalar)*max_row_length*dim_local);
    cudaMallocManaged(&fvec, sizeof(Scalar)*N);
    cudaMallocManaged(&result, sizeof(Scalar)*dim_local);
    cudaMallocManaged(&indices, sizeof(int)*max_row_length*dim_local);

    gpu_ax<<<dim_local,max_row_length>>>(data,fvec,result,indices);
    cudaDeviceSyncrhonize();

//    float elapsed = timer.stop();
//    printf("unified memory takes %f ms to complete with max row length %i and dim local %i \n", elapsed,max_row_length,dim_local);
}
template void mult_vec_unified<int>(int* data, int* fvec, int* result, int* indices, int max_row_length, int dim_local);
template void mult_vec_unified<float>(float* data, float* fvec, float* result, int* indices, int max_row_length, int dim_local);
template void mult_vec_unified<double>(double* data, double* fvec, double* restult, int* indices, int max_row_length, int dim_local);


template<typename Scalar>
void mult_vec_zero(Scalar *data, Scalar *fvec, Scalar *result, int *indices, int max_row_length, int dim_local)
{
    Scalar *data_map, *fvec_map, *result_map;
    int *indices_map;

    cudaGetDeviceProperties(&prop,0);
    if(prop.canMapHostMemory)
    {
      cudaSetDeviceFlags(cudaDeviceMapHost);

      cudaHostAlloc(&data, sizeof(Scalar)*max_row_length*dim_local, cudaHostAllocMapped);
      cudaHostAlloc(&fvec, sizeof(Scalar)*N, cudaHostAllocMapped);
      cudaHostAlloc(&result, sizeof(Scalar)*dim_local, cudaHostAllocMapped);
      cudaHostAlloc(&indices, sizeof(int)*max_row_length*dim_local, cudaHostAllocMapped);

      cudaHostGetDevicePointer(&data_map, data, 0);
      cudaHostGetDevicePointer(&fvec_map, fvec, 0);
      cudaHostGetDevicePointer(&result_map, result, 0);
      cudaHostGetDevicePointer(&indices_map, indices, 0);

      gpu_ax<<<dim_local,max_row_length>>>(data_map,fvec_map,result_map,indices_map);

//      float elapsed = timer.stop();
//      printf("zero copy takes %f ms to complete with max row length %i and dim local %i \n", elapsed,max_row_length,dim_local);
    }
}
template void mult_vec_zero<int>(int* data, int* fvec, int* result, int* indices, int max_row_length, int dim_local);
template void mult_vec_zero<float>(float* data, float* fvec, float* result, int* indices, int max_row_length, int dim_local);
template void mult_vec_zero<double>(double* data, double* fvec, double* restult, int* indices, int max_row_length, int dim_local);

