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
__global__ void kernel(type *vector, type *data, int *indices, type *result, int dim, int max_row_length)
{
    int idx = threadIdx.x + blockDim.x*blockIdx.x;
    if (idx<dim)
    {
        int col;
        type svalue = 0, value;
        for (int i = 0; i < max_row_length; i++)
        {
            value = data[i*dim_local + idx];
            col = indices[i*dim_local + idx];
            svalue += value*vector[col];
        }
        result[idx] = svalue;
    }
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
   *num_blocks = ceil((double)dim / 1024);
   *num_threads = ceil(((double)dim / *num_blocks) / 32) * 32;
}

//=============================================================================
///////////////////////////////////////////////////////////////////////////////
///                             KERNEL TIMING                               ///
///////////////////////////////////////////////////////////////////////////////                       
//=============================================================================
template<typename type>
float invoke_gpu_time(type *vector ,type *data, int * indices, type *result, int dim, int max_row_length, int runs)
{
    Timer timer;
    float elapsed_time = 0.0;
    int num_threads, num_blocks;
    generate_config(&num_threads, &num_blocks, dim);

    //=================================//
    timer.start();
    for (int i = 0; i < runs; i++)
    {
        kernel<<<num_blocks, num_threads>>>(vector, data, indices, result, dim, max_row_length);
    }
    cudaDeviceSynchronize();
    elapsed_time = timer.stop()*1.0e3;
    //=================================//
            
    return elapsed_time / runs;
}
template float invoke_gpu_time<float>(float *vector, float *data, int * indices, float *result, int dim, int max_row_length, int runs);
template float invoke_gpu_time<double>(double *vector, double *data, int * indices, double *result, int dim, int max_row_length, int runs);


//=============================================================================
///////////////////////////////////////////////////////////////////////////////
///                             KERNEL                                      ///
///////////////////////////////////////////////////////////////////////////////                       
//=============================================================================
template<typename type>
void invoke_gpu_overall(type *vector, type *data, int * indices, type *result, int dim, int max_row_length)
{
    int num_threads, num_blocks;
    generate_config(&num_threads, &num_blocks, dim);

    type *placehold = NULL;
    cudaMallocManaged((void **)&placehold, sizeof(type)*num_blocks);

    kernel<<<num_blocks, num_threads, sizeof(double)*num_threads >>>(vector, data, indices, result, dim, max_row_length);
    
    cudaDeviceSynchronize();
}
template void invoke_gpu_overall<float>(float *vector, float *data, int * indices, float *result, int dim, int max_row_length);
template void invoke_gpu_overall<double>(double *vector, double *data, int * indices, double *result, int dim, int max_row_length);



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
