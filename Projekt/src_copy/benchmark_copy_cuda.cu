#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <string.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cmath>
#include <string>
#include "include/timer.hpp"
#define RESET "\e[0m"
#define BLUE "\e[34;1m"
#define CYAN "\e[36;1m"
#define GREY "\e[30;1m"
#define MAGENTA "\e[35;1m"

template <typename Scalar>
void cleanup(Scalar *pointer, int method);

//KERNEL!!!
template<typename type>
__global__ void gpu_copy(type *one, type *two, int dim_local)
{
    int idx = threadIdx.x + blockIdx.x*blockDim.x;
    if (idx < dim_local)
    {
        two[idx] = one[idx];
    }
    
}

//CHANGE!!!!
template<typename type>
void performance(float time_ku, float time_ou, float time_kz, float time_oz, int runs, type schalter, int meth, int ver_first, int ver_second, int mem_option, int dim_local)
{
    using std::string;
    string first, second, method;
    string memop = "";
    if (mem_option == 0) { memop = "(Unified Memory)"; }
    else { memop = "(Zero Copy)"; }
    
    if (meth == 0)
    {
        method = "Unified Memory vs Zero Copy";
        first = "Unified Memory";
        second = "Zero Copy";
    }
    if (meth == 1)
    {
        method = "Kernel vs Kernel";
        first = "Kernel Version" + memop + ": " ;
        first += ver_first;
        second = "Kernel Version" + memop + ": ";
        second += ver_second;

    }
    if (meth == 2)
    {
        method = "Kernel vs CPU";
        first = "Kernel Version" + memop + ": ";
        first += ver_first;
        second = "CPU";
    }
    
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop,0);

    //CHANGE BELOW HERE!!!
    //===#ELEMENTS IN THE VECTOR===================================//
    unsigned long long int elements = dim_local;

    //==='DISK STORAGE~============================================//
    unsigned long long int storage = sizeof(type)*dim_local;
    
    //===#FLOP=====================================================//
    int num_threads = 1024;
    if (dim_local<1024)
    {
        int n = dim_local-1;
        n |= n >> 1;
        n |= n >> 2;
        n |= n >> 4;
        n |= n >> 8;
        n |= n >> 16;
        n |= n >> 16;
        n |= n >> 16;
        num_threads = n + 1;
    }
    int num_blocks = ceil((double)dim_local / 1024);
    
    unsigned long long int flop = elements;         //MULT INTO SHARED MEMORY
    for(int i=2;i<num_threads;i=i*2)
    {
        elements += i*num_blocks;                   //REDUCE KERNEL
    }
    elements += num_blocks;                         //PLACEHOLDER REDUCE

    //==#BYTES=====================================================//           
    int bytes = (elements+1)*sizeof(type);


    printf(GREY "===============================================\n");
    printf(MAGENTA "                PERFORMANCE\n");
    printf("           %s\n", method.c_str());
    printf("        DIM = %i ~~ %i Iterations\n", dim_local, runs);
    printf("            %.2fGB/2GB DRAM used\n", storage / 1.0e9);
    printf(GREY "===============================================\n");
    printf("-----------------------------------------------\n");
    printf(CYAN "                    %s\n", first.c_str());
    printf(GREY "-----------------------------------------------\n");
    printf(CYAN "Kernel Runtime:\t\t\t%f(ms)\n",time_ku);
    printf("Overall Runtime:\t\t%f(ms)\n",time_ou*1.0e3);
    printf("Bandwith(th. Peak):\t\t%.2f(14.9)(GB/s)\n", bytes / (time_ku*1.0e6));
    printf("Flops(th. Peak):\t\t%.6f(326)(GFLOPS/s)\n", flop  / (time_ku*1.0e6));
    printf(GREY "-----------------------------------------------\n");
    printf("-----------------------------------------------\n");
    printf(BLUE "                     %s\n", second.c_str());
    printf(GREY "-----------------------------------------------\n");
    printf(BLUE "Kernel Runtime:\t\t\t%f(ms)\n",time_kz);
    printf("Overall Runtime:\t\t%f(ms)\n",time_oz*1.0e3);
    printf("Bandwith(th. Peak):\t\t%.2f(14.9)(GB/s)\n", bytes / (time_kz*1.0e6));
    printf("Flops(th. Peak):\t\t%.6f(326)(GFLOPS/s)\n", flop  / (time_kz*1.0e6));
    printf(GREY "-----------------------------------------------\n" RESET);



    
}
template void performance<int>(float time_ku, float time_ou, float time_kz, float time_oz, int runs, int schalter, int meth, int ver_first, int ver_second, int mem_option, int dim_local);
template void performance<float>(float time_ku, float time_ou, float time_kz, float time_oz, int runs, float schalter, int meth, int ver_first, int ver_second, int mem_option, int dim_local);
template void performance<double>(float time_ku, float time_ou, float time_kz, float time_oz, int runs, double schalter, int meth, int ver_first, int ver_second, int mem_option, int dim_local);


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

//=============================================================================
//                          ALLOCATION
//                    0=UNIFIED ~~ 1=ZERO COPY
//=============================================================================
template<typename Scalar>
void allocation(Scalar **vector_in, Scalar **vector_out, int dim_local, int mem_option)
{
    switch (mem_option)
    {
    case(0):
        cudaMallocManaged((void **)vector_in, sizeof(Scalar)*dim_local);
        cudaMallocManaged((void **)vector_out, sizeof(Scalar)*dim_local);
        break;
    case(1):
        cudaSetDeviceFlags(cudaDeviceMapHost);
        cudaHostAlloc((void **)vector_in, sizeof(Scalar)*dim_local, cudaHostAllocMapped);
        cudaHostAlloc((void **)vector_out, sizeof(Scalar)*dim_local, cudaHostAllocMapped);
        break;
    }   
}
template void allocation<int>(int **vector_in, int **vector_out, int dim_local, int mem_option);
template void allocation<float>(float **vector_in, float **vector_out, int dim_local, int mem_option);
template void allocation<double>(double **vector_in, double **vector_out, int dim_local, int mem_option);


//=============================================================================
///////////////////////////////////////////////////////////////////////////////
//                          KERNEL
///////////////////////////////////////////////////////////////////////////////
//=============================================================================
template<typename Scalar>
float gpu_copy_time(Scalar *one, Scalar * two, Scalar *result, int dim_local, int runs, int version, int mem_option)
{
    Timer timer;
    float elapsed_time = 0.0;

    int num_blocks = ceil((double)dim_local / 1024);
    int num_threads = ceil(((double)dim_local / num_blocks) / 32) * 32;

    switch (version)
    {
    case(0) :               //kernel_standart
        if (mem_option == 0)
        {
            //=================================//
            timer.start();
            for (int i = 0; i < runs; i++)
            {
                gpu_copy<<<num_blocks, num_threads>>>(one, two, dim_local);
            }
            cudaDeviceSynchronize();
            elapsed_time = timer.stop()*1.0e3;
            //=================================//
        }
        else if (mem_option == 1)
        {
            Scalar *d_one, *d_two;
            cudaHostGetDevicePointer((void **)&d_one, (void *)one, 0);
            cudaHostGetDevicePointer((void **)&d_two, (void *)two, 0);

            //=================================//
            timer.start();
            for (int i = 0; i < runs; i++)
            {
                gpu_copy << <num_blocks, num_threads, sizeof(double)*num_threads >> >(d_one, d_two, dim_local);
            }
            cudaDeviceSynchronize();
            elapsed_time = timer.stop()*1.0e3;
            //=================================//
            
        }
        break;
    }
    return elapsed_time / runs;
}
template float gpu_copy_time<int>(int *one, int * two, int dim_local, int runs, int version, int mem_option);
template float gpu_copy_time<float>(float *one, float * two, int dim_local, int runs, int version, int mem_option);
template float gpu_copy_time<double>(double *one, double * two, int dim_local, int runs, int version, int mem_option);



//GENERATING KERNEL TIME UNIFIED MEMORY
template<typename Scalar>
void gpu_copy_overall(Scalar *one, Scalar * two, int dim_local, int version, int mem_option)
{
    int num_blocks = ceil((double)dim_local / 1024);
    int num_threads = ceil(((double)dim_local / num_blocks) / 32) * 32;

    switch (version)
    {
    case(0) :               //kernel_standart
        if (mem_option == 0)
        {
            gpu_copy << <num_blocks, num_threads >> >(one, two, dim_local);
            cudaDeviceSynchronize();
        }
        else if (mem_option == 1)
        {
            cudaHostGetDevicePointer((void **)&d_one, (void *)one, 0);
            cudaHostGetDevicePointer((void **)&d_two, (void *)two, 0);

            gpu_copy << <num_blocks, num_threads >> >(d_one, d_two, dim_local);
            cudaDeviceSynchronize();
        }
    }
}
template void gpu_copy_overall<int>(int *one, int * two, int dim_local, int version, int mem_option);
template void gpu_copy_overall<float>(float *one, float * two, int dim_local, int version, int mem_option);
template void gpu_copy_overall<double>(double *one, double * two, int dim_local, int version, int mem_option);


//=============================================================================
//                              CLEANUP FUNCTIONS
//=============================================================================
template <typename Scalar>
void cleanup(Scalar *one, Scalar *two, int method)
{
    switch(method)
    {
        case(0):
            cudaFree(one);
            cudaFree(two);
            break;
        case(1):
            cudaFreeHost(one);
            cudaFreeHost(two);
            break;
        case(2):
            delete[] one;
            delete[] two;
            break;
    }
}
template void cleanup<int>(int *one, int *two, int method);
template void cleanup<float>(float *one, float *two, int method);
template void cleanup<double>(double *one, double *two, int method);








//====Ich war nicht mutig genug es zu loeschen :D===/
/*

template <typename Scalar>
void cleanupgpu(Scalar *data)
{
cudaFreeHost(data);
}
template void cleanupgpu<int>(int *data);
template void cleanupgpu<float>(float *data);
template void cleanupgpu<double>(double *data);

//ALLOCATE MEMORY FUNCTION FOR UNIFIED MEMORY for DistEllpack
template<typename Scalar>
void alloc_unifiedD(Scalar **data, int **indices, int dim_local, int dim_local)
{
cudaMallocManaged((void **)data, sizeof(Scalar)*dim_local*dim_local);
cudaMallocManaged((void **)indices, sizeof(int)*dim_local*dim_local);
}
template void alloc_unifiedD<int>(int **data, int **indices, int dim_local, int dim_local);
template void alloc_unifiedD<float>(float **data, int **indices, int dim_local, int dim_local);
template void alloc_unifiedD<double>(double **data, int **indices, int dim_local, int dim_local);

// ALLOCATE MEMORY FUNCTION FOR UNIFIED MEMORY FOR SLICEDVECTOR
template<typename Scalar>
void alloc_unifiedV(Scalar **fvec, int dim_fvec)
{
cudaMallocManaged((void **)fvec, sizeof(Scalar)*dim_fvec);
}
template void alloc_unifiedV<int>(int **fvec, int dim_fvec);
template void alloc_unifiedV<float>(float **fvec, int dim_fvec);
template void alloc_unifiedV<double>(double **fvec, int dim_fvec);
*/