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
__global__ void gpu_axpy(type *scalar, type *vector_x, type *vector_y, type *result, int dim)
{
    idx = threadIdx.x + blockIdx.x*blockDim.x;
    type scale = scalar[0];
    if (idx < dim)
    {
        result[idx] = scale*vector_x[idx] + vector_y[idx];
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
void generate_config(int *num_threads, int *num_blocks)
{

}

//=============================================================================
///////////////////////////////////////////////////////////////////////////////
///                             KERNEL TIMING                               ///
///////////////////////////////////////////////////////////////////////////////                       
//=============================================================================
template<typename type>
float invoke_gpu_time(type scalar, type *vector_x, type *vector_y, type *result, int dim, int runs)
{
    Timer timer;
    float elapsed_time = 0.0;
    int num_threads, num_blocks;
    generate_config(&num_threads, &num_blocks);

    //=================================//
    timer.start();
    for (int i = 0; i < runs; i++)
    {
        gpu_axpy<<<num_blocks, num_threads>>>(scalar, vector_x, vector_y, result, dim);
    }
    cudaDeviceSynchronize();
    elapsed_time = timer.stop()*1.0e3;
    //=================================//
            
    return elapsed_time / runs;
}
template float invoke_gpu_time<float>(float scalar, float *vector_x, float *vector_y, float *result, int dim, int runs);
template float invoke_gpu_time<double>(double scalar, double *vector_x, double *vector_y, double *result, int dim, int runs);


//=============================================================================
///////////////////////////////////////////////////////////////////////////////
///                             KERNEL                                      ///
///////////////////////////////////////////////////////////////////////////////                       
//=============================================================================
template<typename type>
void invoke_gpu_overall(type scalar, type *vector_x, type *vector_y, type *result, int dim)
{
    int num_threads, num_blocks;
    generate_config(&num_threads, &num_blocks);

    gpu_axpy<<<num_blocks, num_threads>>>();
    cudaDeviceSynchronize();
}
template void invoke_gpu_overall<float>(float scalar, float *vector_x, float *vector_y, float *result, int dim);
template void invoke_gpu_overall<double>(double scalar, double *vector_x, double *vector_y, double *result, int dim);



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
