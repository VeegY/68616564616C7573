#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <string.h>
#include <cuda.h>
#include <cuda_runtime.h>
using namespace std;

__global__ void fib_gpu(int *fib)
{
    int temp=fib[1];
    fib[1]+=fib[0];
    fib[0]=temp;
}

extern  void Fibonacci(int* h_fib)
{


     
    
     int *d_fib;

     if(cudaSuccess != cudaMalloc(&d_fib,sizeof(int)*2))
     {
          cout << "allocate error" << endl;
     }
     
     if(cudaSuccess != cudaMemcpy(d_fib, h_fib, sizeof(int)*2, cudaMemcpyHostToDevice))
     {
          cout << "copy error" << endl;
     }

     fib_gpu<<<1, 1>>>(d_fib);
        
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

  cudaDeviceSynchronize();
}
