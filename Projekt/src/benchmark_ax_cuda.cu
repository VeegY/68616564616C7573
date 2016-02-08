#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <string.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cmath>
#include "include/timer.hpp"

template <typename Scalar>
void cleanup(Scalar *data, Scalar *fvec, Scalar *result, int *indices, int method);

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


//CALCULATING MEMORY BANDWITH
template<typename type>
void performance(int max_row_length, int dim_local, float time_ku, float time_ou, float time_kz, float time_oz, int runs, type schalter)
{
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop,0);
    //===#FLOP=====================================================//
    unsigned int flop = 7 * dim_local - 2 - 2 * (floor(pow(dim_local, (1.0 / 3.0)))) - 2 * (floor(pow(dim_local, (2.0 / 3.0)))); //Number of Elements
    flop *= 2;

    //===Immer selbstständig updaten wenn sich der Kernel ändert===//           
    int bRead = 0, bWrite = 0, bytes;
    bRead += max_row_length*dim_local*sizeof(type); //data-Array
    bRead += max_row_length*dim_local*sizeof(int);  //indices-Array
    bRead += max_row_length*dim_local*sizeof(type); //fvec-Array
    bWrite += dim_local*sizeof(type);               //result-Array
    bytes = bRead + bWrite;

    printf("===============================================\n");
    printf("                PERFORMANCE\n");
    printf("          DIM = %i ~~ %i Iterations\n", dim_local, runs);
    printf("===============================================\n");
    printf("-----------------------------------------------\n");
    printf("                UNIFIED_MERMORY\n");
    printf("-----------------------------------------------\n");
    printf("Kernel Runtime:\t\t\t%f(ms)\n",time_ku);
    printf("Overall Runtime:\t\t%f(ms)\n",time_ou*1.0e3);
    printf("Bandwith(th. Peak):\t\t%.2f(14.9)(GB/s)\n", bytes / ((time_ku*1.0e-3)*1.0e9));
    printf("Flops(th. Peak):\t\t%.2f(326)(GFLOPS/s)\n", flop  / ((time_ku*1.0e-3)*1.0e9));
    printf("-----------------------------------------------\n");
    printf("-----------------------------------------------\n");
    printf("                ZERO_COPY\n");
    printf("-----------------------------------------------\n");
    printf("Kernel Runtime:\t\t\t%f(ms)\n",time_kz);
    printf("Overall Runtime:\t\t%f(ms)\n",time_oz*1.0e3);
    printf("Bandwith(th. Peak):\t\t%.2f(14.9)(GB/s)\n", bytes / ((time_kz*1.0e-3)*1.0e9));
    printf("Flops(th. Peak):\t\t%.2f(326)(GFLOPS/s)\n", flop  / ((time_kz*1.0e-3)*1.0e9));
    printf("-----------------------------------------------\n");



    
}
template void performance<int>(int max_row_length, int dim_local, float time_ku, float time_ou, float time_kz, float time_oz, int runs, int schalter);
template void performance<float>(int max_row_length, int dim_local, float time_ku, float time_ou, float time_kz, float time_oz, int runs, float schalter);
template void performance<double>(int max_row_length, int dim_local, float time_ku, float time_ou, float time_kz, float time_oz, int runs, double schalter);


//PROPERTIES OF TEGRA K1
void print_p()
{

    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop,0);
    
    printf("==============================\nDevice name: %s\n------------------------------\n", prop.name);
    printf("Memory Clock Rate (KHz): %d\n",prop.memoryClockRate);
    printf("Memory Bus Width (bits): %d\n",prop.memoryBusWidth);
    printf("Peak Memory Bandwidth (GB/s): %f\n==============================\n",2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/1.0e6);
    //printf("Max Threads Per Block: %i\n", prop.maxThreadsPerBlock);
    //printf("Max Grid Size: %ix%ix%i\n", prop.maxGridSize[0],prop.maxGridSize[1],prop.maxGridSize[2]);

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


//ALLOCATE MEMORY FUNCTION FOR UNIFIED MEMORY for DistEllpack
template<typename Scalar>
void alloc_unifiedD(Scalar **data, int **indices, int max_row_length, int dim_local)
{
    cudaMallocManaged((void **)data, sizeof(Scalar)*dim_local*max_row_length);
    cudaMallocManaged((void **)indices, sizeof(int)*dim_local*max_row_length);
}
template void alloc_unifiedD<int>(int **data, int **indices, int max_row_length, int dim_local);
template void alloc_unifiedD<float>(float **data, int **indices, int max_row_length, int dim_local);
template void alloc_unifiedD<double>(double **data, int **indices, int max_row_length, int dim_local);

// ALLOCATE MEMORY FUNCTION FOR UNIFIED MEMORY FOR SLICEDVECTOR
template<typename Scalar>
void alloc_unifiedV(Scalar **fvec, int dim_fvec)
{
    cudaMallocManaged((void **)fvec, sizeof(Scalar)*dim_fvec);
}
template void alloc_unifiedV<int>(int **fvec, int dim_fvec);
template void alloc_unifiedV<float>(float **fvec, int dim_fvec);
template void alloc_unifiedV<double>(double **fvec, int dim_fvec);


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

//=============================================================================
//                          UNIFIED KERNEL FUNCTIONS
//=============================================================================

//GENERATING KERNEL TIME UNIFIED MEMORY
template<typename Scalar>
float mult_vec_unified_time(Scalar *data, Scalar *fvec, Scalar *result, int *indices, int max_row_length, int dim_local, int dim_fvec, int runs)
{
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    
    int num_blocks = ceil((double)dim_local/1024);
    int num_threads = ceil(((double)dim_local/num_blocks)/32)*32;
    
    cudaEventRecord(start);
    for (int i = 0; i < runs; i++)
    {
        gpu_ax<<<num_blocks,num_threads>>>(data,fvec,result,indices,max_row_length, dim_local);
        
    }
    cudaDeviceSynchronize();
    
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float elapsedTime = 0;
    cudaEventElapsedTime(&elapsedTime, start, stop);

    return (elapsedTime / (float)runs);
}
template float mult_vec_unified_time<int>(int* data, int* fvec, int* result, int* indices, int max_row_length, int dim_local,int dim_fvec, int runs);
template float mult_vec_unified_time<float>(float* data, float* fvec, float* result, int* indices, int max_row_length, int dim_local, int dim_fvec, int runs);
template float mult_vec_unified_time<double>(double* data, double* fvec, double* restult, int* indices, int max_row_length, int dim_local, int dim_fvec, int runs);


//GENERATING KERNEL TIME UNIFIED MEMORY
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


//=============================================================================
//                              ZERO KERNEL FUNCTIONS
//=============================================================================

//KERNE CALL WITH ZERO COPY (NEED TO CALL ALLOC_ZERO BEFORE)
template<typename Scalar>
float mult_vec_zero_time(Scalar *data, Scalar *fvec, Scalar *result, int *indices, int max_row_length, int dim_local, int dim_fvec, int runs)
{
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    Scalar *d_data, *d_fvec, *d_result;
    int *d_indices;

    cudaHostGetDevicePointer((void **)&d_data,(void *)data, 0);
    cudaHostGetDevicePointer((void **)&d_fvec, (void *)fvec, 0);
    cudaHostGetDevicePointer((void **)&d_result, (void *)result, 0);
    cudaHostGetDevicePointer((void **)&d_indices, (void *)indices, 0);

    int num_blocks = ceil((double)dim_local/1024);
    int num_threads = ceil(((double)dim_local/num_blocks)/32)*32;

    cudaEventRecord(start);

    for (int i=0;i<runs;i++)
    {
        gpu_ax<<<num_blocks,num_threads>>>(d_data, d_fvec, d_result, d_indices, max_row_length, dim_local);
   
    }
    cudaDeviceSynchronize();
    
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float elapsedTime = 0;
    cudaEventElapsedTime(&elapsedTime, start, stop);

    cleanup(d_data, d_fvec, d_result, d_indices, 0);
    return (elapsedTime/(float)runs);
}
template float mult_vec_zero_time<int>(int* data, int* fvec, int* result, int* indices, int max_row_length, int dim_local, int  dim_fvec, int runs);
template float mult_vec_zero_time<float>(float* data, float* fvec, float* result, int* indices, int max_row_length, int dim_local, int dim_fvec, int runs);
template float mult_vec_zero_time<double>(double* data, double* fvec, double* restult, int* indices, int max_row_length, int dim_local, int dim_fvec, int runs);


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
    cleanup(d_data, d_fvec, d_result, d_indices, 0);
}
template void mult_vec_zero<int>(int* data, int* fvec, int* result, int* indices, int max_row_length, int dim_local, int  dim_fvec);
template void mult_vec_zero<float>(float* data, float* fvec, float* result, int* indices, int max_row_length, int dim_local, int dim_fvec);
template void mult_vec_zero<double>(double* data, double* fvec, double* restult, int* indices, int max_row_length, int dim_local, int dim_fvec);


//=============================================================================
//                              CLEANUP FUNCTIONS
//=============================================================================
template <typename Scalar>
void cleanup(Scalar *data, Scalar *fvec, Scalar *result, int *indices, int method)
{
    switch(method)
    {
        case(0):
            cudaFree(data);
            cudaFree(fvec);
            cudaFree(result);
            cudaFree(indices);
            break;
        case(1):
            cudaFreeHost(data);
            cudaFreeHost(fvec);
            cudaFreeHost(result);
            cudaFreeHost(indices);
            break;
        case(2):
            delete[] data;
            delete[] fvec;
            delete[] result;
            delete[] indices;
            break;
    }
}
template void cleanup<int>(int *data, int *fvec, int *result, int *indices, int method);
template void cleanup<float>(float *data, float *fvec, float *result, int *indices, int method);
template void cleanup<double>(double *data, double *fvec, double *result, int *indices, int method);


template <typename Scalar>
void cleanupgpu(Scalar *data)
{   
    cudaFreeHost(data);
}
template void cleanupgpu<int>(int *data);
template void cleanupgpu<float>(float *data);
template void cleanupgpu<double>(double *data);
