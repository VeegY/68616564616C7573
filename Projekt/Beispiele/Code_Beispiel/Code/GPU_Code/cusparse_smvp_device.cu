// compile nvcc cusparse_smvp_device.cu -arch=sm_20 -Xcompiler=-fopenmp -lcusparse -o test

// ################################################################################################
//				  Header: cusparse_smvp_device
// ------------------------------------Doxygen-Dokumentation---------------------------------------
///  \file cusparse_smvp_device.cu
///  \brief
///  sparse Matrix (im CSR Format) Vektor Produkt, ausgeführt auf dem Device mit Hilfe von Lib: cuSparse
// ------------------------------------------------------------------------------------------------
// ################################################################################################

#include "helper_functions.hpp"
#include "csr_helper_functions.hpp"
#include <iostream>
#include "cusparse.h"
#define dim 4

void print(double A[], int  m, int n, int lda)
{
	for(int i(0);i<n;++i)
	{
		for (int j(0);j<m;++j)
		{
			std::cout << A[(j*lda)+i] << "  ";
		}
		std::cout << std::endl;
	}
}
void loeschen(cusparseHandle_t handle, cusparseMatDescr_t descr, int *nnzPerRowA, double *A, int *csrRowPtrA, int *csrColIndA, double *csrValA)
{
cusparseDestroy(handle);
cusparseDestroyMatDescr(descr); 
cudaFree(A);
cudaFree(csrValA);
cudaFree(csrRowPtrA);
cudaFree(csrColIndA);
cudaFree(nnzPerRowA);
}

int main()
{
 
  cusparseStatus_t status;
  const int m_dim_matrix(dim),n_dim_matrix(dim),lda(dim);
  double  h_matrix[m_dim_matrix*n_dim_matrix]= {1,0,2,0,3,4,5,0,0,6,7,0,0,0,8,9}; 
  std::cout << "Die Matrix A_h: " << std::endl;
  print(h_matrix,m_dim_matrix,n_dim_matrix,lda);
  
  
  double *d_matrix;
  int *nnzPerRow; 
  int nonzero = 0;
  
 
  
  cudaMalloc(&d_matrix,sizeof(double)*m_dim_matrix*n_dim_matrix);
  cudaMalloc(&nnzPerRow,sizeof(int)*n_dim_matrix);
  
  cudaMemcpy(d_matrix,&h_matrix,sizeof(double)*n_dim_matrix*m_dim_matrix, cudaMemcpyHostToDevice);
  
  
  
  cusparseHandle_t handle=NULL;
  cusparseMatDescr_t descr=NULL;
  cusparseCreate(&handle);
  cusparseCreateMatDescr(&descr); 
  cusparseSetMatType(descr,CUSPARSE_MATRIX_TYPE_GENERAL); // Beschreibt Allgemeine Matrix (z.B. könnte die Matrix symmetrisch sein und so wir ein anderer algorithmus verwendet)
  cusparseSetMatIndexBase(descr,CUSPARSE_INDEX_BASE_ZERO);// Setzt ob der Erste Matrix eintrag a11 ist oder a00
  
  status =  cusparseDnnz(handle, CUSPARSE_DIRECTION_ROW, m_dim_matrix, n_dim_matrix, descr, d_matrix, lda, nnzPerRow, &nonzero);
 
   if (status != CUSPARSE_STATUS_SUCCESS)
   {
     std::cout << "Error "<<std::endl;
     loeschen(handle, descr, nnzPerRow, d_matrix, NULL, NULL, NULL);
     return 0;
   }
  std::cout << "Anzahl der nicht null Elemente  " << nonzero << std::endl;
    
    
  double *csrValA;
  int *csrRowPtrA;
  int *csrColIndA;
  
  
  cudaMalloc(&csrValA,sizeof(float)*nonzero);
  cudaMalloc(&csrRowPtrA,sizeof(int)*(n_dim_matrix+1));
  cudaMalloc(&csrColIndA,sizeof(int)*nonzero);
  
status = cusparseDdense2csr(handle, m_dim_matrix, n_dim_matrix, descr, d_matrix, lda, nnzPerRow, csrValA, csrRowPtrA, csrColIndA);
  if (status != CUSPARSE_STATUS_SUCCESS)
   {
     std::cout << "Error "<<std::endl;
     loeschen(handle, descr, nnzPerRow, d_matrix, csrRowPtrA, csrColIndA, csrValA);
     return 0;
   }
   
cudaDeviceSynchronize();
  
  
    
float alpha(1);
float beta(1);
float h_x[m_dim_matrix]={1,2,3,4};
float h_y[m_dim_matrix];
float *d_x;
float *d_y;



  
cudaMalloc(&d_x,sizeof(float)*m_dim_matrix);
cudaMalloc(&d_y,sizeof(float)*m_dim_matrix);
cudaMemcpy(d_x,&h_x,sizeof(float)*m_dim_matrix, cudaMemcpyHostToDevice);
cudaMemcpy(d_y,&h_y,sizeof(float)*m_dim_matrix, cudaMemcpyHostToDevice);
  
  

 status = cusparseScsrmv((cusparseHandle_t)handle,
			 (cusparseOperation_t)CUSPARSE_OPERATION_TRANSPOSE,
			 (int)m_dim_matrix,
			 (int)n_dim_matrix,
			 (int)nonzero,
			 (const float *) &alpha,
			 (const cusparseMatDescr_t)descr,
			 (const float *)csrValA,
			 (const int *)csrRowPtrA,
			 (const int *)csrColIndA,
			 (const float *)d_x,
			 (const float *)&beta,
			 (float *)d_y);

  if (status != CUSPARSE_STATUS_SUCCESS)
   {
     std::cout << "Error "<<std::endl;
     loeschen(handle, descr, nnzPerRow, d_matrix, csrRowPtrA, csrColIndA, csrValA);
     cudaFree(d_y);
     cudaFree(d_x);
     return 0;
   }

     
  cudaMemcpy(&h_y,d_y,sizeof(float)*m_dim_matrix, cudaMemcpyDeviceToHost);
  
std::cout << std::endl;
std::cout << "y: " << std::endl;
//print(X,n,m,ldc);
print_float(h_y,m_dim_matrix);
std::cout << std::endl;


  
  
  
  // Befreie reservierten Speicher, handle, beschreibung der Matrix
  loeschen(handle, descr, nnzPerRow, d_matrix, csrRowPtrA, csrColIndA, csrValA);
    
}
