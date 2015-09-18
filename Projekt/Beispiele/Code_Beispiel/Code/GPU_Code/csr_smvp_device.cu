// ################################################################################################
//					 FD - Code Project
// ################################################################################################
//				  Header: csr_smvp_device
// ------------------------------------Doxygen-Dokumentation---------------------------------------
///  \file smvp_device.cu
///  \brief
///  Sparse (im CSR Format) Matrix-Vektor Produkt, ausgeführt auf dem Device
// -------------------------------------------------------------------------------------------------
// Kompilieren mit dem folgenden Compilerbefehl: !! WIP !!
//   nvcc csr_smvp_device.cu -O3 -arch=sm_20 -Xcompiler=-fopenmp,-Wall -o prog
//   nvcc csr_smvp_device.cu -arch=sm_20 -o prog
// ------------------------------------------------------------------------------------------------
// Verwendete Header:
#	include <iostream>
#	include <cstdlib>
#	include <ctime>
#	include <cuda_runtime.h>

#	include "helper_functions.hpp"
#	include "csr_helper_functions.hpp"
#	define dim 4
using namespace std;
// #################################################################################################

// ========================================DOKUMENTATION============================================
//				       csr_smvp_oneperrow_sf
///\brief	Sparse Matrix-Vektor Produkt, ein Warp pro Reihe, Matrix im CSR Format
// -------------------------------------------------------------------------------------------------
///@param	value     	1D Array der Größe "nicht-Null-Einträge" der Matrix auf der das CSR Format beruht     
///@param	column		1D Array der selben Größe wie "value", beinhalted die Spalten Indizes der nicht-Null-Einträgen
///@param	pointer   	1D Array der Größe "dim+1" beinhalted den Zeiger auf das erste Element jeder Zeile in "value"
///@param 	vector    	1D Array der Größe dim, dient zur Berechnung des Produktes
///@param 	solution 	1D Array der Größe dim, beinhalted nach Ausführung das Ergebniss
// =================================================================================================
__global__ void csr_smvp_warpperrow_sf(float* value, int* column, int* pointer, float* vector,float* solution)
{
     int idx=threadIdx.x;
     solution[idx%5]=1;
     if(idx<pointer[blockIdx.x+1])
     {
          value[pointer[blockIdx.x]+idx]*=vector[column[pointer[blockIdx.x]+idx]];
          //atomicAdd(&solution[blockIdx.x],value[pointer[blockIdx.x]+idx]);
          //atomicAdd(&solution[blockIdx.x],1);
          
          
     }

}
// =================================================================================================

// =================================================================================================
//					 main-Funktion
// -------------------------------------------------------------------------------------------------
///  \brief Fuehrt den Test der csr_smvp_device Implementierung durch.
// =================================================================================================
int main()
{

     /*========Variabeln======*/

     float *h_matrix= new float[dim*dim];
     float *h_vector= new float[dim];
     float *h_solution = new float[dim];
     float *d_value;
     float *d_solution;
     float *d_vector;
     
     int *d_column;
     int *d_pointer;
     
     
     sparse_mat_fill(h_matrix,2,dim);
     dense_vec_fill(h_vector,dim);
     
        
     int nz= nz_values(h_matrix,dim);
     cout << nz << endl << "------------------------------------------------------------------------------" << endl;
     int *h_column = new int[nz];
     int *h_pointer= new int[dim+1];
     
     float *h_value = new float[nz];
 
     print_float(h_vector,dim);
     csr_format_cpu(h_matrix,h_value,h_column,h_pointer,dim);
     print_mat(h_matrix,dim);
     print_float(h_value,nz);
     print_int(h_column,nz);
     print_int(h_pointer,dim+1);
     

     

    
     if(cudaSuccess != cudaMalloc(&d_value, sizeof(float)*nz))
          {
               cout << "allocate error" << endl;
          }
     
     if(cudaSuccess != cudaMalloc(&d_column, sizeof(int)*nz))
          {
               cout << "allocate error" << endl;
               cudaFree(d_value);
          }     
     
     if(cudaSuccess != cudaMalloc(&d_pointer, sizeof(int)*(dim+1)))
          {
               cout << "allocate error" << endl;
               cudaFree(d_value);
               cudaFree(d_column);
          } 
          
      if(cudaSuccess != cudaMalloc(&d_solution, sizeof(float)*dim))
          {
               cout << "allocate error" << endl;
               cudaFree(d_value);
               cudaFree(d_column);
               cudaFree(d_pointer);
          } 
      if(cudaSuccess != cudaMalloc(&d_vector, sizeof(float)*dim))
          {
               cout << "allocate error" << endl;
               cudaFree(d_value);
               cudaFree(d_column);
               cudaFree(d_pointer);
               cudaFree(d_solution);
          } 
           
         
     
     
     /*========Memcpy=========*/    
     if(cudaSuccess != cudaMemcpy(d_value, h_value, sizeof(float)*nz, cudaMemcpyHostToDevice))
          {
               cout << "failed to copy" << endl;
               cudaFree(d_value);
               cudaFree(d_column);
               cudaFree(d_pointer);
               cudaFree(d_solution);
               cudaFree(d_vector);
          }
          
     if(cudaSuccess != cudaMemcpy(d_column, h_column, sizeof(int)*nz, cudaMemcpyHostToDevice))
          {
               cout << "failed to copy" << endl;
               cudaFree(d_value);
               cudaFree(d_column);
               cudaFree(d_pointer);
               cudaFree(d_solution);
               cudaFree(d_vector);
          }          
     
     if(cudaSuccess != cudaMemcpy(d_pointer, h_pointer, sizeof(int)*(dim+1), cudaMemcpyHostToDevice))
          {
               cout << "failed to copy" << endl;
               cudaFree(d_value);
               cudaFree(d_column);
               cudaFree(d_pointer);
               cudaFree(d_solution);
               cudaFree(d_vector);
          }     
     
     if(cudaSuccess != cudaMemcpy(d_solution, h_solution, sizeof(float)*dim, cudaMemcpyHostToDevice))
          {
               cout << "failed to copy" << endl;
               cudaFree(d_value);
               cudaFree(d_column);
               cudaFree(d_pointer);
               cudaFree(d_solution);
               cudaFree(d_vector);
          } 
     
     if(cudaSuccess != cudaMemcpy(d_vector, h_vector, sizeof(float)*dim, cudaMemcpyHostToDevice))
          {
               cout << "failed to copy" << endl;
               cudaFree(d_value);
               cudaFree(d_column);
               cudaFree(d_pointer);
               cudaFree(d_solution);
               cudaFree(d_vector);
          }
     /*=======Kernel==========*/
     csr_smvp_warpperrow_sf<<<dim,32>>>(d_value,d_column,d_pointer,d_vector,d_solution);
     if(cudaSuccess != cudaGetLastError())
          {
               cout << "kernel launch failed" << endl;
               cudaFree(d_value);
               cudaFree(d_column);
               cudaFree(d_pointer);
               cudaFree(d_solution);
               cudaFree(d_vector);   
          }
          
          
     /*=======Memcpy==========*/
     if(cudaSuccess != cudaMemcpy(h_solution, d_solution, sizeof(float)*dim, cudaMemcpyDeviceToHost))
          {
               cout << "failed to copy - device to host" << endl;
               cudaFree(d_value);
               cudaFree(d_column);
               cudaFree(d_pointer);
               cudaFree(d_solution);
               cudaFree(d_vector);
          }     
     
    print_float(h_solution,dim);
    
     
     /*========Cleanup========*/
     cudaFree(d_value);
     cudaFree(d_column);
     cudaFree(d_pointer);
     cudaFree(d_solution);
     cudaFree(d_vector);
     
     delete[] h_matrix;
     delete[] h_vector;
     delete[] h_solution;
     
     return 0;
     
}
// =================================================================================================
