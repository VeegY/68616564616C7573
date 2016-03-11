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
__global__ void kernel()
{

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
float invoke_gpu_time(int runs)
{
    Timer timer;
    float elapsed_time = 0.0;
    int num_threads, num_blocks;
    generate_config(&num_threads, &num_blocks);

    //=================================//
    timer.start();
    for (int i = 0; i < runs; i++)
    {
        kernel<<<num_blocks, num_threads>>>();
    }
    cudaDeviceSynchronize();
    elapsed_time = timer.stop()*1.0e3;
    //=================================//
            
    return elapsed_time / runs;
}
template float invoke_gpu_time<float>(int runs);
template float invoke_gpu_time<double>(int runs);


//=============================================================================
///////////////////////////////////////////////////////////////////////////////
///                             KERNEL                                      ///
///////////////////////////////////////////////////////////////////////////////                       
//=============================================================================
template<typename type>
void invoke_gpu_overall()
{
    int num_threads, num_blocks;
    generate_config(&num_threads, &num_blocks);

    kernel<<<num_blocks, num_threads>>>();
    cudaDeviceSynchronize();
}
template void invoke_gpu_overall<float>();
template void invoke_gpu_overall<double>();



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
