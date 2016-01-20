#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <string.h>
#include <cuda.h>
#include <cuda_runtime.h>

//KERNEL
template<typename type>
__global__ void  gpu_ax_opt(type* data, type* fvec, type* result, int* indices, int dim_local)
{
	/*
	Kernel Idee:
	-extern shared memory für data oder data und indices -testen!
	-values im shared array speicher, d/i
	-berrechnen über alle Werte (data = 0 auslassen)? hinderlich für memory coalescence?)
	-danach reduce Teil über max_row_length auf gleichen shared memory
	-ein Block hat die Größe: ceil: (shared_memory/32) keine halbe Reihen!
	*/
	/*
		Block Größe Detail
		shared memory: 32kB
		Datentyp:
		 foat(4),int(4),double(8)
		Menge:
		 Type*max_row_length*dim_local
		+int*max_row_length*dim_local
		wähle große Blöcke, keine halben Reihen und vielfaches von 32!
     */
}


//KERNEL
template<typename type>
__global__ void  gpu_ax(type* data, type* fvec, type* result, int* indices, int max_row_length, int dim_local)
{

	int idx = 0;
	bool zero = true;
    type value = 0;
    while(zero)
    {
		if (data[idx*dim_local+threadIdx.x] == 0)
		{
			zero = false;
		}
		else
		{
			value += data[idx*dim_local+threadIdx.x] * fvec[indices[idx*dim_local+threadIdx.x]];
			idx++;
		}
      
    }
    result[threadIdx.x]=value;
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


//KERNEL CALL WITH UNIFIED MEMORY (NEED TO CALL ALLOC_UNIFIED BEFORE)
template<typename Scalar>
void mult_vec_unified(Scalar *data, Scalar *fvec, Scalar *result, int *indices, int max_row_length, int dim_local, int dim_fvec)
{

    gpu_ax<<<1,max_row_length>>>(data,fvec,result,indices,max_row_length, dim_local);
    cudaDeviceSynchronize();

}
template void mult_vec_unified<int>(int* data, int* fvec, int* result, int* indices, int max_row_length, int dim_local,int dim_fvec);
template void mult_vec_unified<float>(float* data, float* fvec, float* result, int* indices, int max_row_length, int dim_local, int dim_fvec);
template void mult_vec_unified<double>(double* data, double* fvec, double* restult, int* indices, int max_row_length, int dim_local, int dim_fvec);


//KERNE CALL WITH ZERO COPY (NEED TO CALL ALLOC_ZERO BEFORE)
template<typename Scalar>
void mult_vec_zero(Scalar *data, Scalar *fvec, Scalar *result, int *indices, int max_row_length, int dim_local, int dim_fvec)
{
    Scalar *d_data, *d_fvec, *d_result;
    int *d_indices;

    cudaHostGetDevicePointer((void **)&d_data,(void *)data, 0);
    cudaHostGetDevicePointer((void **)&d_fvec, (void *)fvec, 0);
    cudaHostGetDevicePointer((void **)&d_result, (void *)result, 0);
    cudaHostGetDevicePointer((void **)&d_indices, (void *)indices, 0);

    gpu_ax<<<1,dim_local>>>(d_data, d_fvec, d_result, d_indices, max_row_length, dim_local);
    cudaDeviceSynchronize();
}
template void mult_vec_zero<int>(int* data, int* fvec, int* result, int* indices, int max_row_length, int dim_local, int  dim_fvec);
template void mult_vec_zero<float>(float* data, float* fvec, float* result, int* indices, int max_row_length, int dim_local, int dim_fvec);
template void mult_vec_zero<double>(double* data, double* fvec, double* restult, int* indices, int max_row_length, int dim_local, int dim_fvec);

template <typename Scalar>
void cleanup(Scalar *data, Scalar *fvec, Scalar *result, int *indices)
{
	cudaFreeHost(data);
	cudaFreeHost(fvec);
	cudaFreeHost(result);
	cudaFreeHost(indices);
}
template void cleanup<int>(int *data, int *fvec, int *result, int *indices);
template void cleanup<float>(float *data, float *fvec, float *result, int *indices);
template void cleanup<double>(double *data, double *fvec, double *result, int *indices);
