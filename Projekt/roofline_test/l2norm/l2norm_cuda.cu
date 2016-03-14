#include <iomanip>
#include <stdio.h>
#include <string.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cmath>
#include "include/timer.hpp"

//=============================================================================
///////////////////////////////////////////////////////////////////////////////
///                             KERNEL                                      ///
///////////////////////////////////////////////////////////////////////////////
//=============================================================================
template<typename type>
__global__ void kernel(type *vector, type *placehold, int dim_local)
{
    extern __shared__ double array[];
    type* shar = (type*)array;
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int sidx = threadIdx.x;
    type value = (type)0;
    if (idx < dim_local)
    {
        value = vector[idx];
        value *= value;
    }
    shar[sidx] = value;
    __syncthreads();

    //reduce kernel
    for (int offset = blockDim.x / 2; offset >0; offset >>= 1)
    {
        if (sidx < offset)
        {
            shar[sidx] += shar[sidx + offset];
        }
        __syncthreads();
    }

    if (sidx == 0)
    {
        placehold[blockIdx.x] = shar[0];
    }
}

template<typename type>
__global__ void resultreduce(type *result, type *placehold, int num_blocks)
{
    type value = (type)0;
    for (int i = 0; i < num_blocks; i++)
    {
        value += placehold[i];
    }
    result[0] = sqrt(value);
}
//=============================================================================
///////////////////////////////////////////////////////////////////////////////
///                             ALLOCATION                                  ///
///////////////////////////////////////////////////////////////////////////////
//=============================================================================
template<typename type>
void allocation(type **data, size_t size)
{    
    cudaMallocManaged((void **)data, sizeof(type)*size);
}
template void allocation<int>(int **data, size_t size);
template void allocation<float>(float **data, size_t size);
template void allocation<double>(double **data, size_t size);


//=============================================================================
///////////////////////////////////////////////////////////////////////////////
///                             KERNEL CONFIG                               ///
///////////////////////////////////////////////////////////////////////////////                       
//=============================================================================
void generate_config(int *num_threads, int *num_blocks, int dim)
{
    
    *num_threads = 1024;
    if (dim<1024)
    {
        int n = dim - 1;
        n |= n >> 1;
        n |= n >> 2;
        n |= n >> 4;
        n |= n >> 8;
        n |= n >> 16;
        n |= n >> 16;
        n |= n >> 16;
        *num_threads = n + 1;
    }
    *num_blocks = ceil((double)dim / 1024);
}

//=============================================================================
///////////////////////////////////////////////////////////////////////////////
///                             KERNEL TIMING                               ///
///////////////////////////////////////////////////////////////////////////////                       
//=============================================================================
template<typename type>
float invoke_gpu_time(type *vector, type *result, int dim, int runs)
{
    Timer timer;
    float elapsed_time = 0.0;
    int num_threads, num_blocks;
    generate_config(&num_threads, &num_blocks, dim);

    type *placehold = NULL;
    cudaMallocManaged((void **)&placehold, sizeof(type)*num_blocks);

    //=================================//
    timer.start();
    for (int i = 0; i < runs; i++)
    {
        kernel<<<num_blocks, num_threads, sizeof(double)*num_threads>>>(vector,placehold,dim);
        resultreduce<<<1, 1>>>(result, placehold, num_blocks);
    }
    cudaDeviceSynchronize();
    elapsed_time = timer.stop()*1.0e3;
    //=================================//
            
    return elapsed_time / runs;
}
template float invoke_gpu_time<float>(float *vector, float *result, int dim, int runs);
template float invoke_gpu_time<double>(double *vector, double *result, int dim, int runs);


//=============================================================================
///////////////////////////////////////////////////////////////////////////////
///                             KERNEL                                      ///
///////////////////////////////////////////////////////////////////////////////                       
//=============================================================================
template<typename type>
void invoke_gpu_overall(type *vector, type *result, int dim)
{
    int num_threads, num_blocks;
    generate_config(&num_threads, &num_blocks, dim);

    type *placehold = NULL;
    cudaMallocManaged((void **)&placehold, sizeof(type)*num_blocks);

    kernel<<<num_blocks, num_threads, sizeof(double)*num_threads >>>(vector, placehold, dim);
    resultreduce << <1, 1 >> >(result, placehold, num_blocks);
    
    cudaDeviceSynchronize();
}
template void invoke_gpu_overall<float>(float *vector, float *result, int dim);
template void invoke_gpu_overall<double>(double *vector, double *result, int dim);



//=============================================================================
///////////////////////////////////////////////////////////////////////////////
///                             CLEANUP                                     ///
///////////////////////////////////////////////////////////////////////////////                       
//=============================================================================
template <typename type>
void cleanup(type *data)
{
   cudaFree(data);
}
template void cleanup<int>(int *data);
template void cleanup<float>(float *data);
template void cleanup<double>(double *data);
