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
__global__ void gpu_maxnorm(type *vector, type *placehold, int dim_local, int numblocks)
{
    extern __shared__ double array[];
    type* shar = (type*)array;
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int sidx = threadIdx.x;
    type value = 0, compare_one, compare_two;
    if (idx < dim_local)
    {
        
        value = vector[idx];
        if (value < 0)
        {
            value = -value;
        }
    }
    shar[sidx] = value;
    __syncthreads();

    //reduce kernel
    for (int offset = blockDim.x / 2; offset >0; offset >>= 1)
    {
        if (sidx < offset)
        {
            compare_one = shar[idx];
            compare_two = shar[idx + offset];
                if (compare_two > compare_one)
                {
                    shar[idx] = compare_two;
                }
        }
        __syncthreads();
    }

    if (sidx == 0)
    {
        placehold[blockIdx.x] = shar[0];
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
void allocation(Scalar **vector, Scalar **result, int dim_local, int mem_option)
{
    switch (mem_option)
    {
    case(0):
        cudaMallocManaged((void **)vector, sizeof(Scalar)*dim_local);
        cudaMallocManaged((void **)result, sizeof(Scalar));
        break;
    case(1):
        cudaSetDeviceFlags(cudaDeviceMapHost);
        cudaHostAlloc((void **)vector, sizeof(Scalar)*dim_local, cudaHostAllocMapped);
        cudaHostAlloc((void **)result, sizeof(Scalar), cudaHostAllocMapped);
        break;
    }   
}
template void allocation<int>(int **vector, int **result, int dim_local, int mem_option);
template void allocation<float>(float **vector, float **result, int dim_local, int mem_option);
template void allocation<double>(double **vector, double **result, int dim_local, int mem_option);


//=============================================================================
//                          KERNEL
//=============================================================================
template<typename Scalar>
float gpu_maxnorm_time(Scalar *vector, Scalar *result, int dim_local, int runs, int version, int mem_option)
{
    Timer timer;
    float elapsed_time = 0.0;

    Scalar *placehold = NULL;
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
    //printf("%f -- %f", num_threads, num_blocks);

    switch (version)
    {
    case(0) :               //kernel_standart
        if (mem_option == 0)
        {
            cudaMallocManaged((void **)&placehold, sizeof(Scalar)*num_blocks);
            //=================================//
            timer.start();
            for (int i = 0; i < runs; i++)
            {
                gpu_maxnorm<<<num_blocks, num_threads, sizeof(double)*num_threads >>>(vector, placehold, dim_local, num_blocks);
            }
            cudaDeviceSynchronize();
            Scalar value = placehold[0];
            for(int i=1;i<num_blocks;i++)
            {
                value += placehold[i];
            }
            result[0] = value;
            elapsed_time = timer.stop()*1.0e3;
            //=================================//
        }
        else if (mem_option == 1)
        {
            cudaHostAlloc((void **)&placehold, sizeof(Scalar)*num_blocks, cudaHostAllocMapped);
            Scalar *d_vector, *d_placehold;

            cudaHostGetDevicePointer((void **)&d_vector, (void *)vector, 0);
            cudaHostGetDevicePointer((void **)&d_placehold, (void *)placehold, 0);

            //=================================//
            timer.start();
            for (int i = 0; i < runs; i++)
            {
                gpu_maxnorm << <num_blocks, num_threads, sizeof(double)*num_threads >> >(d_vector, d_placehold, dim_local, num_blocks);
            }
            cudaDeviceSynchronize();
            Scalar value = placehold[0];
            for(int i=1;i<num_blocks;i++)
            {
                if (value < placehold[i])
                {
                    value = placehold[i];
                }
            }
            result[0] = sqrt(value);
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
template float gpu_maxnorm_time<int>(int *vector, int *result, int dim_local, int runs, int version, int mem_option);
template float gpu_maxnorm_time<float>(float *vector, float *result, int dim_local, int runs, int version, int mem_option);
template float gpu_maxnorm_time<double>(double *vector, double *result, int dim_local, int runs, int version, int mem_option);



//GENERATING KERNEL TIME UNIFIED MEMORY
template<typename Scalar>
void gpu_maxnorm_overall(Scalar *vector, Scalar *result, int dim_local, int version, int mem_option)
{
    Scalar *placehold = NULL;

    int num_threads = 1024;
    if (dim_local<1024)
    {
        int n = dim_local - 1;
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
    switch (version)
    {
    case(0) :               //kernel_standart
        if(mem_option == 0)
        {
            cudaMallocManaged((void **)&placehold, sizeof(Scalar)*num_blocks);
            gpu_maxnorm <<<num_blocks, num_threads, sizeof(double)*num_threads>>>(vector, placehold, dim_local, num_blocks);
            cudaDeviceSynchronize();
            Scalar value = placehold[0];
            for(int i=1;i<num_blocks;i++)
            {
                value += placehold[i];
            }
            result[0] = sqrt(value);
        }
        else if(mem_option == 1)
        {
            cudaHostAlloc((void **)&placehold, sizeof(Scalar)*num_blocks, cudaHostAllocMapped);
            Scalar *d_vector, *d_placehold;

            cudaHostGetDevicePointer((void **)&d_vector, (void *)vector, 0);
            cudaHostGetDevicePointer((void **)&d_placehold, (void *)placehold, 0);
            
            gpu_maxnorm<<<num_blocks, num_threads, sizeof(double)*num_threads>>>(d_vector, d_placehold, dim_local, num_blocks);
            cudaDeviceSynchronize();
            Scalar value = placehold[0];
            for(int i=1;i<num_blocks;i++)
            {
                value += placehold[i];
            }
            result[0] = sqrt(value);
        }
        break;

    case(1) :               //kernel_shared(NUR ALS BEISPIEL)
            //gpu_ax_shared << <num_blocks, num_threads >> >(d_data, d_fvec, d_result, d_indices, dim_local, dim_local);
        break;

    case(2) :               //kernel_advanced(NUR ALS BEISPIEL)
           //gpu_ax_advanced << <num_blocks, num_threads >> >(d_data, d_fvec, d_result, d_indices, dim_local, dim_local);
        break;
    }
}
template void gpu_maxnorm_overall<int>(int *vector, int *result, int dim_local, int version, int mem_option);
template void gpu_maxnorm_overall<float>(float *vector, float *result, int dim_local, int version, int mem_option);
template void gpu_maxnorm_overall<double>(double *vector, double *result, int dim_local, int version, int mem_option);


//=============================================================================
//                              CLEANUP FUNCTIONS
//=============================================================================
template <typename Scalar>
void cleanup(Scalar *vector, Scalar *result, int method)
{
    switch(method)
    {
        case(0):
            cudaFree(vector);
            cudaFree(result);
            break;
        case(1):
            cudaFreeHost(vector);
            cudaFreeHost(result);
            break;
        case(2):
            delete[] vector;
            delete[] result;
            break;
    }
}
template void cleanup<int>(int *vector, int *result, int method);
template void cleanup<float>(float *vector, float *result, int method);
template void cleanup<double>(double *vector, double *result, int method);
