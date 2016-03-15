#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <string.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cmath>
#include <string>
//#include "include/timer.hpp"


template <typename Scalar>
void cleanup(Scalar *data, Scalar *fvec, Scalar *result, int *indices, int method);

//KERNEL
template<typename mtype, typename vtype, typename rtype>
__global__ void  gpu_ax(mtype* data, const vtype* fvec, rtype* result, size_t* indices, size_t max_row_length, size_t dim_local)
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


template<typename type>
void performance(int max_row_length, int dim_local, float time_ku, float time_ou, float time_kz, float time_oz, int runs, type schalter, int meth, int ver_first, int ver_second, int mem_option)
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

    //===#ELEMENTS IN THE MATRIX===================================//
    unsigned long long int elements = 7 * dim_local - 2 - 2 * (floor(pow(dim_local, (1.0 / 3.0)))) - 2 * (floor(pow(dim_local, (2.0 / 3.0))));

    //==='DISK STORAGE~============================================//
    unsigned long long int storage = sizeof(type)*(2 * dim_local + dim_local*max_row_length) + sizeof(int)*dim_local*max_row_length;

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
template void performance<int>(int max_row_length, int dim_local, float time_ku, float time_ou, float time_kz, float time_oz, int runs, int schalter, int meth, int ver_first, int ver_second, int mem_option);
template void performance<float>(int max_row_length, int dim_local, float time_ku, float time_ou, float time_kz, float time_oz, int runs, float schalter, int meth, int ver_first, int ver_second, int mem_option);
template void performance<double>(int max_row_length, int dim_local, float time_ku, float time_ou, float time_kz, float time_oz, int runs, double schalter, int meth, int ver_first, int ver_second, int mem_option);


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
void allocation(Scalar **data, Scalar **fvec, Scalar **result, int **indices, int max_row_length, int dim_local,int dim_fvec, int mem_option)
{
    switch (mem_option)
    {
    case(0):
        cudaMallocManaged((void **)data, sizeof(Scalar)*dim_local*max_row_length);
        cudaMallocManaged((void **)fvec, sizeof(Scalar)*dim_fvec);
        cudaMallocManaged((void **)result, sizeof(Scalar)*dim_local);
        cudaMallocManaged((void **)indices, sizeof(int)*dim_local*max_row_length);
        break;
    case(1):
        cudaSetDeviceFlags(cudaDeviceMapHost);
        cudaHostAlloc((void **)data, sizeof(Scalar)*max_row_length*dim_local, cudaHostAllocMapped);
        cudaHostAlloc((void **)fvec, sizeof(Scalar)*dim_fvec, cudaHostAllocMapped);
        cudaHostAlloc((void **)result, sizeof(Scalar)*dim_local, cudaHostAllocMapped);
        cudaHostAlloc((void **)indices, sizeof(int)*max_row_length*dim_local, cudaHostAllocMapped);
        break;
    }
}
template void allocation<int>(int **data, int **fvec, int **result, int **indices, int max_row_length, int dim_local, int dim_fvec, int mem_option);
template void allocation<float>(float **data, float **fvec, float **result, int **indices, int max_row_length, int dim_local, int dim_fvec, int mem_option);
template void allocation<double>(double **data, double **fvec, double **result, int **indices, int max_row_length, int dim_local, int dim_fvec, int mem_option);


//=============================================================================
//                          KERNEL
//=============================================================================
template<typename Scalar>
float gpu_ax_time(Scalar *data, Scalar *fvec, Scalar *result, int *indices, int max_row_length, int dim_local, int dim_fvec, int runs, int version, int mem_option)
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
                gpu_ax << <num_blocks, num_threads >> >(data, fvec, result, indices, max_row_length, dim_local);
            }
            cudaDeviceSynchronize();
            elapsed_time = timer.stop()*1.0e3;
            //=================================//
        }
        else if (mem_option == 1)
        {
            Scalar *d_data, *d_fvec, *d_result;
            int *d_indices;

            cudaHostGetDevicePointer((void **)&d_data, (void *)data, 0);
            cudaHostGetDevicePointer((void **)&d_fvec, (void *)fvec, 0);
            cudaHostGetDevicePointer((void **)&d_result, (void *)result, 0);
            cudaHostGetDevicePointer((void **)&d_indices, (void *)indices, 0);

            //=================================//
            timer.start();
            for (int i = 0; i < runs; i++)
            {
                gpu_ax << <num_blocks, num_threads >> >(d_data, d_fvec, d_result, d_indices, max_row_length, dim_local);
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
template float gpu_ax_time<int>(int* data, int* fvec, int* result, int* indices, int max_row_length, int dim_local, int dim_fvec, int runs, int version, int mem_option);
template float gpu_ax_time<float>(float* data, float* fvec, float* result, int* indices, int max_row_length, int dim_local, int dim_fvec, int runs, int version, int mem_option);
template float gpu_ax_time<double>(double* data, double* fvec, double* restult, int* indices, int max_row_length, int dim_local, int dim_fvec, int runs, int version, int mem_option);



//=============================================================================
//                          KERNEL für DistEllpackKlasse
//=============================================================================
template<typename mtype, template vtype, template rtype>
void gpu_ax_(mtype* data, const vtype* fvec, rtype* result, size_t *indices, size_t max_row_length,
		  size_t dim_local)
{
    int num_blocks = ceil((double)dim_local / 1024);
    int num_threads = ceil(((double)dim_local / num_blocks) / 32) * 32;

            //=================================//
                gpu_ax << <num_blocks, num_threads >> >(data, fvec, result, indices, max_row_length, dim_local);
            cudaDeviceSynchronize();
            //=================================//
}
template void gpu_ax_<float, float, float>(mtype*, const vtype*, rtype*, size_t*, size_t, size_t);
template void gpu_ax_<double, double, double>(mtype*, const vtype*, rtype*, size_t*, size_t, size_t);



//GENERATING KERNEL TIME UNIFIED MEMORY
template<typename Scalar>
void gpu_ax_overall(Scalar *data, Scalar *fvec, Scalar *result, int *indices, int max_row_length, int dim_local, int dim_fvec, int version, int mem_option)
{
    int num_blocks = ceil((double)dim_local / 1024);
    int num_threads = ceil(((double)dim_local / num_blocks) / 32) * 32;

    switch (version)
    {
    case(0) :               //kernel_standart
        if(mem_option == 0)
        {
            gpu_ax <<<num_blocks, num_threads >>>(data, fvec, result, indices, max_row_length, dim_local);
        }
        else if(mem_option == 1)
        {
            Scalar *d_data, *d_fvec, *d_result;
            int *d_indices;

            cudaHostGetDevicePointer((void **)&d_data, (void *)data, 0);
            cudaHostGetDevicePointer((void **)&d_fvec, (void *)fvec, 0);
            cudaHostGetDevicePointer((void **)&d_result, (void *)result, 0);
            cudaHostGetDevicePointer((void **)&d_indices, (void *)indices, 0);

            gpu_ax <<<num_blocks, num_threads >>>(d_data, d_fvec, d_result, d_indices, max_row_length, dim_local);
        }
        cudaDeviceSynchronize();
        break;

    case(1) :               //kernel_shared(NUR ALS BEISPIEL)
            //gpu_ax_shared << <num_blocks, num_threads >> >(d_data, d_fvec, d_result, d_indices, max_row_length, dim_local);
        break;

    case(2) :               //kernel_advanced(NUR ALS BEISPIEL)
           //gpu_ax_advanced << <num_blocks, num_threads >> >(d_data, d_fvec, d_result, d_indices, max_row_length, dim_local);
        break;
    }
}
template void gpu_ax_overall<int>(int* data, int* fvec, int* result, int* indices, int max_row_length, int dim_local,int dim_fvec, int version, int mem_option);
template void gpu_ax_overall<float>(float* data, float* fvec, float* result, int* indices, int max_row_length, int dim_local, int dim_fvec, int version, int mem_option);
template void gpu_ax_overall<double>(double* data, double* fvec, double* restult, int* indices, int max_row_length, int dim_local, int dim_fvec, int version, int mem_option);


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








//====Ich war nicht mutig genug es zu loeschen :D===/

template <typename Scalar>
void cleanupgpu(Scalar *data)
{
cudaFree(data);
}
template void cleanupgpu<int>(int *data);
template void cleanupgpu<float>(float *data);
template void cleanupgpu<double>(double *data);


//ALLOCATE MEMORY FUNCTION FOR UNIFIED MEMORY for DistEllpack
template<typename Scalar>
void alloc_unifiedD(Scalar **data, size_t **indices, int max_row_length, int dim_local)
{
cudaMallocManaged((void **)data, sizeof(Scalar)*dim_local*max_row_length);
cudaMallocManaged((void **)indices, sizeof(size_t)*dim_local*max_row_length);
}
template void alloc_unifiedD<int>(int **data, size_t **indices, int max_row_length, int dim_local);
template void alloc_unifiedD<float>(float **data, size_t **indices, int max_row_length, int dim_local);
template void alloc_unifiedD<double>(double **data, size_t **indices, int max_row_length, int dim_local);


// ALLOCATE MEMORY FUNCTION FOR UNIFIED MEMORY FOR SLICEDVECTOR
template<typename Scalar>
void alloc_unifiedV(Scalar **fvec, int dim_fvec)
{
cudaMallocManaged((void **)fvec, sizeof(Scalar)*dim_fvec);
}
template void alloc_unifiedV<int>(int **fvec, int dim_fvec);
template void alloc_unifiedV<float>(float **fvec, int dim_fvec);
template void alloc_unifiedV<double>(double **fvec, int dim_fvec);
