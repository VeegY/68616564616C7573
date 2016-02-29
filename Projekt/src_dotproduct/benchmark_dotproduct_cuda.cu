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
__global__ void gpu_dotproduct(type *one,type *two, type *result, int dim_local)
{
    extern __shared__ type shar[];
    idx = threadIdx.x + blockIdx.x * blockDim.x;
    type value = (type)0;
    if (idx < dim_local)
    {
        value = one[idx];
        value *= two[idx];
    }
    shar[idx] = value;
    __syncthreads();

    //reduce kernel
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
        atomicAdd(&result, shar[idx]));
    }

}

/*//CHANGE!!!!
template<typename type>
void performance(float time_ku, float time_ou, float time_kz, float time_oz, int runs, type schalter, int meth, int ver_first, int ver_second, int mem_option)
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
    //===#ELEMENTS IN THE MATRIX===================================//
    unsigned long long int elements = 7 * dim_local - 2 - 2 * (floor(pow(dim_local, (1.0 / 3.0)))) - 2 * (floor(pow(dim_local, (2.0 / 3.0))));

    //==='DISK STORAGE~============================================//
    unsigned long long int storage = sizeof(type)*(2 * dim_local + dim_local*dim_local) + sizeof(int)*dim_local*dim_local;
    
        //===#FLOP=====================================================//
    unsigned long long int flop = 2 * elements;

    //==#BYTES=====================================================//           
    int bytes = elements*(sizeof(type) + sizeof(int)) + 2*(sizeof(type)*dim_local);// Elements(Data+Indices) + Fvec Read und Result Write
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
template void performance<int>(float time_ku, float time_ou, float time_kz, float time_oz, int runs, int schalter, int meth, int ver_first, int ver_second, int mem_option);
template void performance<float>(float time_ku, float time_ou, float time_kz, float time_oz, int runs, float schalter, int meth, int ver_first, int ver_second, int mem_option);
template void performance<double>(float time_ku, float time_ou, float time_kz, float time_oz, int runs, double schalter, int meth, int ver_first, int ver_second, int mem_option);
*/

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
void allocation(Scalar **vecone, Scalar **vectwo, Scalar **result, int dim_local, int mem_option)
{
    switch (mem_option)
    {
    case(0):
        cudaMallocManaged((void **)vecone, sizeof(Scalar)*dim_local);
        cudaMallocManaged((void **)vectwo, sizeof(Scalar)*dim_local);
        cudaMallocManaged((void **)result, sizeof(Scalar));
        break;
    case(1):
        cudaSetDeviceFlags(cudaDeviceMapHost);
        cudaHostAlloc((void **)vecone, sizeof(Scalar)*dim_local, cudaHostAllocMapped);
        cudaHostAlloc((void **)vectwo, sizeof(Scalar)*dim_local, cudaHostAllocMapped);
        cudaHostAlloc((void **)result, sizeof(Scalar), cudaHostAllocMapped);
        break;
    }   
}
template void allocation<int>(int **vecone, int **vectwo, int **result, int dim_local, int mem_option);
template void allocation<float>(float **vecone, float **vectwo, float **result, int dim_local, int mem_option);
template void allocation<double>(double **vecone, double **vectwo, double **result, int dim_local, int mem_option);


//=============================================================================
//                          KERNEL
//=============================================================================
template<typename Scalar>
float gpu_dotproduct_time(Scalar *one, Scalar * two, Scalar *result, int dim_local, int runs, int version, int mem_option)
{
    Timer timer;
    float elapsed_time = 0.0;

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
    printf("%f -- %f", num_threads, num_blocks);

    switch (version)
    {
    case(0) :               //kernel_standart
        if (mem_option == 0)
        {
            //=================================//
            timer.start();
            for (int i = 0; i < runs; i++)
            {
                gpu_dotproduct<<<num_blocks, num_threads, sizeof(Scalar)*dim_local>>>(one, two, result, dim_local);
            }
            cudaDeviceSynchronize();
            elapsed_time = timer.stop()*1.0e3;
            //=================================//
        }
        else if (mem_option == 1)
        {
            Scalar *d_one, *d_two, *d_result;

            cudaHostGetDevicePointer((void **)&d_one, (void *)one, 0);
            cudaHostGetDevicePointer((void **)&d_two, (void *)two, 0);
            cudaHostGetDevicePointer((void **)&d_result, (void *)result, 0);

            //=================================//
            timer.start();
            for (int i = 0; i < runs; i++)
            {
                gpu_dotproduct << <num_blocks, num_threads, sizeof(Scalar)*dim_local >> >(d_one, d_two, d_result, dim_local);
            }
            cudaDeviceSynchronize();
            elapsed_time = timer.stop()*1.0e3;
            //=================================//
            
        }
        break;

    case(1) :   //kernel_shared(NUR ALS BEISPIEL)
                /*//=================================//
                timer.start();
                for (int i = 0; i < runs; i++)
                {
                        kernel_shared<<<>>>
                }
                cudaDeviceSynchronize();
                elapsed_time = timer.stop();
                //=================================//*/
        break;

    case(2) :   //kernel_advanced(NUR ALS BEISPIEL)
                /*//=================================//
                timer.start();
                for (int i = 0; i < runs; i++)
                {
                        kernel_advanced<<<>>>
                }
                cudaDeviceSynchronize();
                elapsed_time = timer.stop();
                //=================================//*/
        break;
    }
    return elapsed_time / runs;
}
template float gpu_dotproduct_time<int>(int *one, int * two, int *result, int dim_local, int runs, int version, int mem_option);
template float gpu_dotproduct_time<float>(float *one, float * two, float *result, int dim_local, int runs, int version, int mem_option);
template float gpu_dotproduct_time<double>(double *one, double * two, double *result, int dim_local, int runs, int version, int mem_option);



//GENERATING KERNEL TIME UNIFIED MEMORY
template<typename Scalar>
void gpu_dotproduct_overall(Scalar *one, Scalar * two, Scalar *result, int dim_local, int version, int mem_option)
{
    int num_threads=1024;
    if(dim_local<1024)
    {
        num_threads = GetNextPowerOf2(dim_local);
    }
    int num_blocks = ceil((double)dim_local / 1024);

    switch (version)
    {
    case(0) :               //kernel_standart
        if(mem_option == 0)
        {
            gpu_dotproduct <<<num_blocks, num_threads, sizeof(Scalar)*dim_local >>>(one, two, result, dim_local);
        }
        else if(mem_option == 1)
        {
            Scalar *d_one, *d_two, *d_result;

            cudaHostGetDevicePointer((void **)&d_one, (void *)one, 0);
            cudaHostGetDevicePointer((void **)&d_two, (void *)two, 0);
            cudaHostGetDevicePointer((void **)&d_result, (void *)result, 0);
            
            gpu_dotproduct << <num_blocks, num_threads, sizeof(Scalar)*dim_local >> >(d_one, d_two, d_result, dim_local);
        }
        cudaDeviceSynchronize();
        break;

    case(1) :               //kernel_shared(NUR ALS BEISPIEL)
            //gpu_ax_shared << <num_blocks, num_threads >> >(d_data, d_fvec, d_result, d_indices, dim_local, dim_local);
        break;

    case(2) :               //kernel_advanced(NUR ALS BEISPIEL)
           //gpu_ax_advanced << <num_blocks, num_threads >> >(d_data, d_fvec, d_result, d_indices, dim_local, dim_local);
        break;
    }
}
template void gpu_dotproduct_overall<int>(int *one, int * two, int *result, int dim_local, int version, int mem_option);
template void gpu_dotproduct_overall<float>(float *one, float * two, float *result, int dim_local, int version, int mem_option);
template void gpu_dotproduct_overall<double>(double *one, double * two, double *result, int dim_local, int version, int mem_option);


//=============================================================================
//                              CLEANUP FUNCTIONS
//=============================================================================
template <typename Scalar>
void cleanup(Scalar *one, Scalar *two, Scalar *result, int method)
{
    switch(method)
    {
        case(0):
            cudaFree(one);
            cudaFree(two);
            cudaFree(result);
            break;
        case(1):
            cudaFreeHost(one);
            cudaFreeHost(two);
            cudaFreeHost(result);
            break;
        case(2):
            delete[] one;
            delete[] two;
            delete[] result;
            break;
    }
}
template void cleanup<int>(int *one, int *two, int *result, int method);
template void cleanup<float>(float *one, float *two, float *result, int method);
template void cleanup<double>(double *one, double *two, double *result, int method);








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