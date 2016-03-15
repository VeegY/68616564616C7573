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
      size_t col;
      rtype svalue = 0;
      mtype value;
      for(int i = 0;i < max_row_length; i++)
      {
        value = data[i*dim_local+idx];
        col = indices[i*dim_local+idx];
        svalue += value*fvec[col];
      }
      result[idx]=svalue;
    }
}



//=============================================================================
//                          KERNEL für DistEllpackKlasse
//=============================================================================
template<typename mtype, typename vtype, typename rtype>
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
template void gpu_ax_<float, float, float>(float*, const float*, float*, size_t*, size_t, size_t);
template void gpu_ax_<double, double, double>(double*, const double*, double*, size_t*, size_t, size_t);




//====Ich war nicht mutig genug es zu loeschen :D===/

template <typename Scalar>
void cleanupgpu(Scalar *data)
{
cudaFree(data);
}
template void cleanupgpu<int>(int *data);
template void cleanupgpu<float>(float *data);
template void cleanupgpu<double>(double *data);
template void cleanupgpu<size_t>(size_t *data);


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
