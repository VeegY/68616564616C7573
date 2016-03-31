#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <string.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cmath>
#include <string>


template<typename Scalar>
void alloc_unified(Scalar **fvec, size_t dim_fvec);

template <typename Scalar>
void cleanupgpu(Scalar *data);


//KERNEL Matrix Vecotr Produkt
template<typename mtype, typename vtype, typename rtype>
__global__ void  gpu_ax(mtype* data, const vtype* fvec, rtype* result, size_t* indices, size_t max_row_length, size_t dim_local)
{

    int idx = threadIdx.x+blockDim.x*blockIdx.x;
    if(idx<dim_local)
    {
      int col;
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



template<typename type> // Dot Produkt Kernel
__global__ void gpu_dot(type *vectorx, type *vectory, type *placehold, int dim_local)
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

template<typename type> // Dot Produkt Kernel mit const
__global__ void gpu_dot(const type *vectorx, const type *vectory, type *placehold, int dim_local)
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




template<typename type> // reduce Kernel für Dot Produkt
__global__ void resultreduce(type *result, type *placehold, int num_blocks)
{
    type value = (type)0;
    for (int i = 0; i < num_blocks; i++)
    {
        value += placehold[i];
    }
    result[0] = value;
}


template<typename type> // Kernel für axpy
__global__ void axpygpu(type *vector_x, type scalar, type *vector_y, type *result, size_t dim)
{
    size_t idx = threadIdx.x + blockIdx.x*blockDim.x;
    if (idx < dim)
    {
        result[idx] = scalar*vector_x[idx] + vector_y[idx];
    }

}



template<typename type> //Kernel für L2-Norm
__global__ void L2_Norm(type *vector, type *placehold, int dim_local)
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

template<typename type> // resultreduce für L2 norm
__global__ void resultreduce_l2(type *result, type *placehold, int num_blocks)
{
    type value = (type)0;
    for (int i = 0; i < num_blocks; i++)
    {
        value += placehold[i];
    }
    result[0] = value;//sqrt(value);
}



template<typename type> // Maxnorm Kernel
__global__ void maxn(type *vector, type *placehold, int dim_local)
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

    //reduce kernel for maxnorm
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



template<typename type>
__global__ void resultreducemax(type *result, type *placehold, int num_blocks)
{
    type value = placehold[0];
    for (int i = 1; i < num_blocks; i++)
    {
        if (value <= placehold[i])
        {
            value = placehold[i];
        }
    }
    result[0] = value;
}



template<typename type> // copy Kernel
__global__ void copygpu(const type *vecin, type scalar, type *vecout, size_t dim)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx < dim)
    {
        vecout[idx] = scalar * vecin[idx];
    }
}


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
//                          Aufruf für DistEllpackKlasse
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





//=============================================================================
//                   Aufruf Dot Produkt (in FillVectorgpu Klasse)            //
//=============================================================================
template<typename type>
void gpu_dot_(type *vecx, type *vecy, size_t dim, type *erg)
{
    int num_threads, num_blocks;
    generate_config(&num_threads, &num_blocks, dim);

    type *placehold = NULL;
//    cudaMallocManaged((void **)&placehold, sizeof(type)*num_blocks);
    cudaMallocManaged(&placehold, sizeof(type)*num_blocks);

    //=================================//
        gpu_dot<<<num_blocks, num_threads, sizeof(type)*num_threads>>>(vecx,vecy,placehold,dim);
        resultreduce<<<1, 1>>>(erg, placehold, num_blocks);

        cudaDeviceSynchronize();

    //=================================//
    cleanupgpu(placehold);
}
template void gpu_dot_<double>(double *vecx, double *vecy, size_t dim, double *erg);
template void gpu_dot_<float>(float *vecx, float *vecy, size_t dim, float *erg);
template void gpu_dot_<int>(int *vecx, int *vecy, size_t dim, int *erg);

//=============================================================================
//                   Aufruf Dot Produkt (in FillVectorgpu Klasse) mit const  //
//=============================================================================
template<typename type>
void gpu_dot_(const type *vecx, const type *vecy, size_t dim, type *erg)
{
    int num_threads, num_blocks;
    generate_config(&num_threads, &num_blocks, dim);

    type *placehold = NULL;
//    cudaMallocManaged((void **)&placehold, sizeof(type)*num_blocks);
    cudaMallocManaged(&placehold, sizeof(type)*num_blocks);

    //=================================//
        gpu_dot<<<num_blocks, num_threads, sizeof(type)*num_threads>>>(vecx,vecy,placehold,dim);
        cudaDeviceSynchronize();//TODO TODELETE
        resultreduce<<<1, 1>>>(erg, placehold, num_blocks);

        cudaDeviceSynchronize();
    //=================================//
    cleanupgpu(placehold);
}
template void gpu_dot_<double>(const double *vecx, const double *vecy, size_t dim, double *erg);
template void gpu_dot_<float>(const float *vecx, const float *vecy, size_t dim, float *erg);
template void gpu_dot_<int>(const int *vecx, const int *vecy, size_t dim, int *erg);

//=============================================================================
//                   Aufruf axpy                                             //
//=============================================================================
template<typename type>
void gpu_axpy(type *vecx, type scalar, type *vecy, size_t dim)
{
    int num_threads, num_blocks;
    generate_config(&num_threads, &num_blocks, dim);


    //=================================//
        axpygpu<<<num_blocks, num_threads, sizeof(type)*num_threads>>>(vecx,scalar,vecy,vecx,dim);
        cudaDeviceSynchronize();

    //=================================//

}
template void gpu_axpy<double>(double *vecx, double scalar, double *vecy, size_t dim);
template void gpu_axpy<float>(float *vecx, float scalar, float *vecy, size_t dim);
template void gpu_axpy<int>(int *vecx, int scalar, int *vecy, size_t dim);


//=============================================================================
//                   Aufruf L2-Norm                                          //
//=============================================================================
template<typename type>
void gpu_l2(type *vec, size_t dim, type *erg)
{
    int num_threads, num_blocks;
    generate_config(&num_threads, &num_blocks, dim);

    type *placehold = NULL;
    cudaMallocManaged((void **)&placehold, sizeof(type)*num_blocks);
//    cudaMallocManaged(&placehold, sizeof(type)*num_blocks);

    //=================================//
        L2_Norm<<<num_blocks, num_threads, sizeof(type)*num_threads>>>(vec,placehold,dim);
        cudaDeviceSynchronize();//TODO TOREMOVE
        resultreduce_l2<<<1, 1>>>(erg, placehold, num_blocks);

        cudaDeviceSynchronize();

    //=================================//

    cleanupgpu(placehold);

}
template void gpu_l2<double>(double *vec, size_t dim, double *erg);
template void gpu_l2<float>(float *vec, size_t dim, float *erg);

//TODO TOADD const


//=============================================================================
//                   Aufruf unendlich-Norm                                   //
//=============================================================================
template<typename type>
void gpumaxnorm(type *vec, size_t dim, type *erg)
{
    int num_threads, num_blocks;
    generate_config(&num_threads, &num_blocks, dim);

    type *placehold = NULL;
    cudaMallocManaged((void **)&placehold, sizeof(type)*num_blocks);
 //   cudaMallocManaged(&placehold, sizeof(type)*num_blocks);

    //=================================//
        maxn<<<num_blocks, num_threads, sizeof(type)*num_threads>>>(vec,placehold,dim);
        resultreducemax<<<1, 1>>>(erg, placehold, num_blocks);

        cudaDeviceSynchronize();

    //=================================//

    cleanupgpu(placehold);

}
template void gpumaxnorm<double>(double *vec, size_t dim, double *erg);
template void gpumaxnorm<float>(float *vec, size_t dim, float *erg);


//=============================================================================
//                          Aufruf für copy Kernel                           //
//=============================================================================
template<typename type>
void copygpu_(const type *vecin, type *vecout, size_t dim)
{
    int num_threads, num_blocks;
    generate_config(&num_threads, &num_blocks, dim);

    type scalar(1);

     //=================================//
        copygpu<<<num_blocks, num_threads, sizeof(type)*num_threads>>>(vecin,scalar,vecout,dim);
        cudaDeviceSynchronize();

    //=================================//
}

template void copygpu_<double>(const double *vecin, double *vecout, size_t dim);
template void copygpu_<float>(const float *vecin, float *vecout, size_t dim);
template void copygpu_<int>(const int *vecin, int *vecout, size_t dim);


//=============================================================================
//                              CLEANUP FUNCTION
//=============================================================================

template <typename Scalar>
void cleanupgpu(Scalar *data)
{
cudaFree(data);
}
template void cleanupgpu<int>(int *data);
template void cleanupgpu<float>(float *data);
template void cleanupgpu<double>(double *data);
template void cleanupgpu<size_t>(size_t *data);

// ALLOCATE MEMORY FUNCTION FOR UNIFIED MEMORY FOR SLICEDVECTOR
template<typename Scalar>
void alloc_unified(Scalar **fvec, size_t dim_fvec)
{
cudaMallocManaged((void **)fvec, sizeof(Scalar)*dim_fvec);
}
template void alloc_unified<int>(int **fvec, size_t dim_fvec);
template void alloc_unified<float>(float **fvec, size_t dim_fvec);
template void alloc_unified<double>(double **fvec, size_t dim_fvec);
template void alloc_unified<size_t>(size_t **fvec, size_t dim_fvec);
