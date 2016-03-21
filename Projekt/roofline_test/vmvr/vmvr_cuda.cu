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
//=================================================================//
template<typename type>
__global__ void kernel_dot(type *vectorx, type *vectory, type *placehold, int dim_local)
{
    extern __shared__ double array[];
    type* shar = (type*)array;
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int sidx = threadIdx.x;
    type value = 0;
    if (idx < dim_local)
    {
        value = vectorx[idx];
        value *= vectory[idx];
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

//=================================================================//
//=================================================================//
template<typename type>
__global__ void kernel_l2norm(type *vector, type *placehold, int dim_local)
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

//=================================================================//
//=================================================================//
template<typename type>
__global__ void resultreduce(type *result, type *placehold, int num_blocks, int loops)
{
    extern __shared__ double array[];
    type* shar = (type*)array;
    int idx = threadIdx.x;
    type value = (type)0;
    
    for (int i = 0; i < loops; i++)
    {
        if (idx+i*1024 < dim_local)
        {
            value += placehold[idx+i*1024];
        }
    }
    shar[idx] = value;
    __syncthreads();
    
    for (int offset = blockDim.x / 2; offset >0; offset >>= 1)
    {
        if (idx < offset)
        {
            shar[idx] += shar[idx + offset];
        }
        __syncthreads();
    }

    if (idx == 0)
    {
        result[0] = sqrt(shar[0]);
    }
}
    
    
    
    /*type value = (double)0;
    for (int i = 0; i < num_blocks; i++)
    {
        value += placehold[i];
    }
    result[0] = sqrt(value);*/
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
float invoke_gpu_time_dotproduct(type *vecx, type *vecy, type *placehold, int dim, int runs)
{
    Timer timer;
    float elapsed_time = 0.0;
    int num_threads, num_blocks;
    generate_config(&num_threads, &num_blocks, dim);
    //=================================//
    timer.start();
    for (int i = 0; i < runs; i++)
    {
        kernel_dot<<<num_blocks, num_threads, sizeof(double)*num_threads>>>(vecx,vecy,placehold,dim);
    }
    cudaDeviceSynchronize();
    elapsed_time = timer.stop()*1.0e3;
    //=================================//
    return elapsed_time / runs;
}
template float invoke_gpu_time_dotproduct<float>(float *vecx, float *vecy, float *placehold, int dim, int runs);
template float invoke_gpu_time_dotproduct<double>(double *vecx, double *vecy, double *placehold, int dim, int runs);


template<typename type>
float invoke_gpu_time_l2norm(type *vec, type *placehold, int dim, int runs)
{
    Timer timer;
    float elapsed_time = 0.0;
    int num_threads, num_blocks;
    generate_config(&num_threads, &num_blocks, dim);

    //=================================//
    timer.start();
    for (int i = 0; i < runs; i++)
    {
        kernel_l2norm<< <num_blocks, num_threads, sizeof(double)*num_threads >> >(vec, placehold, dim);  
    }
    cudaDeviceSynchronize();
    elapsed_time = timer.stop()*1.0e3;
    //=================================//
    return elapsed_time / runs;
}
template float invoke_gpu_time_l2norm<float>(float *vec, float *placehold, int dim, int runs);
template float invoke_gpu_time_l2norm<double>(double *vec, double *placehold, int dim, int runs);

template<typename type>
float invoke_gpu_time_reduce(type *placehold, type *result, int dim, int runs)
{
    Timer timer;
    float elapsed_time = 0.0;
    int num_threads, num_blocks;
    generate_config(&num_threads, &num_blocks, dim);
    int loops = ceil((double)num_blocks / 1024);
    
    //=================================//
    timer.start();
    for (int i = 0; i < runs; i++)
    {
        resultreduce << <1, 1024,sizeof(double)*1024 >> >(result, placehold, num_blocks);
    }
    cudaDeviceSynchronize();
    elapsed_time = timer.stop()*1.0e3;
    //=================================//

    return elapsed_time / runs;
}
template float invoke_gpu_time_reduce<float>(float *placehold, float *result, int dim, int runs);
template float invoke_gpu_time_reduce<double>(double *placehold, double *result, int dim, int runs);


//=============================================================================
///////////////////////////////////////////////////////////////////////////////
///                             KERNEL                                      ///
///////////////////////////////////////////////////////////////////////////////                       
//=============================================================================
template<typename type>
void invoke_gpu_overall(type *vecx, type *vecy, type *result, int dim)
{
    int num_threads, num_blocks;
    generate_config(&num_threads, &num_blocks, dim);

    type *placehold = NULL;
    cudaMallocManaged((void **)&placehold, sizeof(type)*num_blocks);

    kernel_dot<<<num_blocks, num_threads, sizeof(double)*num_threads >>>(vecx, vecy, placehold, dim);
    resultreduce << <1, 1 >> >(result, placehold, num_blocks);
    
    cudaDeviceSynchronize();
}
template void invoke_gpu_overall<float>(float *vecx, float *vecy, float *result, int dim);
template void invoke_gpu_overall<double>(double *vecx, double *vecy, double *result, int dim);



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
