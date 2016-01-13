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

    if(prop.canMapHostMemory)
    {
      cudaSetDeviceFlags(cudaDeviceMapHost);

      std::cout << "IFPART\n";

      int* h_a = NULL;
      int* h_b = NULL;
      int* h_c = NULL;


      cudaHostAlloc((void **)&h_a, sizeof(int)*N, cudaHostAllocMapped);
      cudaHostAlloc((void **)&h_b, sizeof(int)*N, cudaHostAllocMapped);
      cudaHostAlloc((void **)&h_c, sizeof(int)*N, cudaHostAllocMapped);

      int *d_a, *d_b, *d_c;

      h_a[0]=1;
      h_a[1]=2;
      h_b[0]=3;
      h_b[1]=4;
      
      cudaHostGetDevicePointer((void **)&d_a, (void *)h_a,0);
      cudaHostGetDevicePointer((void **)&d_b, (void *)h_b,0);
      cudaHostGetDevicePointer((void **)&d_c, (void *)h_c,0);

      vadd<<<1,N>>>(d_a,d_b,d_c);
      cudaDeviceSynchronize();


     // printf("c[0]=%i\nc[1]=%i\n",c[0],c[1]);
      printf("h_c[0]=%i\nh_c[1]=%i\n",h_c[0],h_c[1]);
     //std::cout << &h_c << " " << h_c << "\n";
    }
    else{std::cout << "CANT MAP HOST MEMORY\n";}

   // printf("c[0]=%i\nc[1]=%i\n",c[0],c[1]);


    return 0;
}

