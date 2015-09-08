/*   CompileInfo: 
     nvcc smpv_device.cu -O3 -arch=sm_20 -Xcompiler=-fopenmp,-Wall -o prog
     nvcc smpv_device.cu -arch=sm_20 -o prog
*/
// ################################################################################################
//				  Header: smvp_device
// ------------------------------------Doxygen-Dokumentation---------------------------------------
///  \file smvp_device.cu
///  \brief
///  sparse Matrix Vektor Produkt, ausgeführt auf dem Device
// ------------------------------------------------------------------------------------------------
// ################################################################################################
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cuda_runtime.h>
#include "helper_functions.hpp"
#define dim 3
#define csr 2
using namespace std;


// ========================================DOKUMENTATION============================================
//					      smvp_oneperrow_sf
///\brief	sparse Matrix Vektor Produkt, ein Thread pro Reihe
// -------------------------------------------------------------------------------------------------
///@param	matrix    	dim x dim Matrix gespeichert als 1D Array der Größe dim*dim
///@param	vector		Vector gespeichert als 1D Array der Größe dim
///@param	solution  	Vektor der Größe dim in dem das Ergebniss gespeichert wird
///@param 	n         	Dimension des Vektor oder Matrix
// =================================================================================================
__global__ void smvp_oneperrow_sf(float* matrix, float* vector, float* solution, uint n)
{
     int idx=threadIdx.x+blockIdx.x*blockDim.x;
     float value =0;
     
     if(idx<n)
     {
          for(int i=0;i<n;i++)
          {
               value += matrix[idx*n+i]*vector[i];
          }          
          solution[idx]=value;
     }
}


// ========================================DOKUMENTATION============================================
//					      smvp_oneperrow_sf
///\brief	sparse Matrix Vektor Produkt, ein Thread pro Element
// -------------------------------------------------------------------------------------------------
///@param	matrix    	dim x dim Matrix gespeichert als 1D Array der Größe dim*dim
///@param	vector		Vector gespeichert als 1D Array der Größe dim
///@param	solution  	Vektor der Größe dim in dem das Ergebniss gespeichert wird
///@param 	n         	Dimension des Vektor oder Matrix
// =================================================================================================
__global__ void smvp_threadperelement_sf(float* matrix, float* vector, float* solution, uint n)
{
     int idx=threadIdx.x+blockIdx.x*blockDim.x;
     if(idx<n*n)
     {
          matrix[idx] *= vector[idx%n];
          atomicAdd(&solution[idx/n],matrix[idx]);
     }     
}


int main()
{

     /*========Variabeln======*/
     
     float *h_matrix = new float[dim*dim];
     float *h_vector= new float[dim];
     float *h_solution = new float[dim];
     float *d_matrix;
     float *d_vector;
     float *d_solution;
     
     sparse_mat_fill(h_matrix,3,dim);
     dense_vec_fill(h_vector,dim);
     
 /*    
     print_float(value,4);
     cout << endl;
     print_int(column,4);
     cout << endl;
     print_int(pointerE,3);
     cout << endl;
 */    
     /*========Malloc=========*/
    
          if(cudaSuccess != cudaMalloc(&d_matrix, sizeof(float)*dim*dim))
          {
               cout << "allocate error" << endl;
          }
     
          if(cudaSuccess != cudaMalloc(&d_vector, sizeof(float)*dim))
          {
               cout << "allocate error" << endl;
               cudaFree(d_matrix);
          }     
     
          if(cudaSuccess != cudaMalloc(&d_solution, sizeof(float)*dim))
          {
               cout << "allocate error" << endl;
               cudaFree(d_matrix);
               cudaFree(d_vector);
          } 
 
     
     
     /*========Memcpy=========*/    
     if(cudaSuccess != cudaMemcpy(d_matrix, h_matrix, sizeof(float)*dim*dim, cudaMemcpyHostToDevice))
          {
               cout << "failed to copy" << endl;
               cudaFree(d_matrix);
               cudaFree(d_vector);
               cudaFree(d_solution);
          }
          
     if(cudaSuccess != cudaMemcpy(d_vector, h_vector, sizeof(float)*dim, cudaMemcpyHostToDevice))
          {
               cout << "failed to copy" << endl;
               cudaFree(d_matrix);
               cudaFree(d_vector);
               cudaFree(d_solution);
          }          
     
     if(cudaSuccess != cudaMemcpy(d_solution, h_solution, sizeof(float)*dim, cudaMemcpyHostToDevice))
          {
               cout << "failed to copy" << endl;
               cudaFree(d_matrix);
               cudaFree(d_vector);
               cudaFree(d_solution);
          }     
     
     
     /*=======Kernel==========*/
     //smvp_oneperrow_sf<<<1,4>>>(d_matrix,d_vector,d_solution,dim);
     smvp_threadperelement_sf<<<dim,dim>>>(d_matrix,d_vector,d_solution,dim);
     
     if(cudaSuccess != cudaGetLastError())
          {
               cout << "kernel launch failed" << endl;
               cudaFree(d_matrix);
               cudaFree(d_vector);
               cudaFree(d_solution);     
          }
          
          
     /*=======Memcpy==========*/
     if(cudaSuccess != cudaMemcpy(h_solution, d_solution, sizeof(float)*dim, cudaMemcpyDeviceToHost))
          {
               cout << "failed to copy" << endl;
               cudaFree(d_matrix);
               cudaFree(d_vector);
               cudaFree(d_solution);
          }     
     
    cout << h_solution[0] << " " << h_solution[1] << " " << h_solution[2] <<endl;
    
     
     /*========Cleanup========*/
     cudaFree(d_matrix);
     cudaFree(d_vector);
     cudaFree(d_solution);
     
     delete[] h_matrix;
     delete[] h_vector;
     delete[] h_solution;
     
     return 0;
     
}
