#include <cuda.h>
#include <cuda_runtime.h>

template <typename Scalar>
__global__ void mult_vec_impl_gpu()
{
    size_t idx=threadIdx.x;
    size_t pos=idx*_dim_local
    Scalar res=0;
    if(idx<_max_row_length)
    {
        
    }
}

template <typename Scalar>
extern void cuda_standart()
{
    Scalar *d_fvec,*d_indices,*d_data,*d_whatelse;
    int size_fvec = N * sizeof(Scalar);
    int size_indices = N * sizeof(int);
    int size_data = N* sizeof(Scalar);
    
    cudaMalloc();
    cudaMalloc();
    cudaMalloc();
    //LOG_ERROR("Cuda: Malloc failed");
    
    cudaMemcpy();
    cudaMemcpy();
    cudaMemcpy();
    //LOG_ERROR("Cuda: Memory copy failed!")
    
    mult_vec_impl_gpu<<<1,1>>>();
}

template <typename Scalar>
extern void cuda_unified()
{
    cudaMallocManaged()
    cudaMallocManaged()
    cudaMallocManaged()
    
    mult_vec_impl_gpu<<<1,1>>>();
        
}


extern void cuda_none()
{

}
