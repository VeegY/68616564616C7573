#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <string.h>
#include <cuda.h>
#include <cuda_runtime.h>
using namespace std;

template <typename Scalar>
__global__ void mult_vec_gpu(Scalar *data, Scalar *fvec, Scalar *result, int *indices, int max_row_length, int dim_local)
{
    
    
}

template <typename Scalar>
extern  void mult_vec_memcpy(Scalar *e_data, Scalar *e_fvec, Scalar *e_result, int *e_indices, int e_max_row_length, int e_dim_local)
{

     Scalar* d_data,d_fvec,d_result;
     int* d_indices;
     
     int size_m = e_max_row_length*e_dim_local;
     cout << "check check\n";
     /*cudaMalloc(&d_data, sizeof(Scalar)*size_m);
     cudaMalloc(&d_fvec, sizeof(Scalar)*e_dim_local); //TODO: dim fvec
     cudaMalloc(&d_result, sizeof(Scalar)*e_dim_local);
     cudaMalloc(&d_indices, sizeof(int)*size_m);
     
     cudaMemcpy(d_data, e_data, sizeof(Scalar)*size_m,cudaMemcpyHostToDevice);
     cudaMemcpy(d_fvec, e_fvec, sizeof(Scalar)*e_dim_local,cudaMemcpyHostToDevice); //TODO: dim fvec 
     cudaMemcpy(d_result, e_result, sizeof(Scalar)*e_dim_local,cudaMemcpyHostToDevice);
     cudaMemcpy(d_indices, e_indices, sizeof(int)*size_m,cudaMemcpyHostToDevice);
     

     fib_gpu<float><<<1, 1>>>(d_data, d_fvec, d_result, d_indices, e_max_row_length, e_dim_local);
        
     if(cudaSuccess != cudaGetLastError())
     {
          cout << "kernel launch failed" << endl;
     }
     
     
     if(cudaSuccess != cudaMemcpy(h_fib, d_fib, sizeof(int)*2, cudaMemcpyDeviceToHost))
     {
          cout << "copy error" << endl;
     }     
//=====================Kernel=============================//

  //GPUtimer gtimer;
  //gtimer.start();

  //double gtime = gtimer.stop();
  //time[1] = gtime;

  cudaDeviceSynchronize();*/
}
