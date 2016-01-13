#include <string.h>
#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <iostream>
#define N 2
__global__ void vadd(int *a, int *b, int *c)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    c[idx]=a[idx]+b[idx];
}

int main()
{
      cudaDeviceProp prop;
      cudaGetDeviceProperties(&prop,0);

      int* h_a = NULL;
      int* h_b = NULL;
      int* h_c = NULL;
      
      h_a[0]=1;
      h_a[1]=2;
      h_b[0]=3;
      h_b[1]=4;

      cudaMallocManaged((void **)&h_a, sizeof(int)*N);
      cudaMallocManaged((void **)&h_b, sizeof(int)*N);
      cudaMallocManaged((void **)&h_c, sizeof(int)*N);

      //h_a[0]=1;
      //h_a[1]=2;
      //h_b[0]=3;
      //h_b[1]=4;

      vadd<<<1,N>>>(h_a,h_b,h_c);
      cudaDeviceSynchronize();

      printf("h_c[0]=%i\nh_c[1]=%i\n",h_c[0],h_c[1]);

    return 0;
}

