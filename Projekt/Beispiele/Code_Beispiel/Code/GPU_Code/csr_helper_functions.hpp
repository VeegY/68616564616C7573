// ################################################################################################
//				  Header: csr_helper_functions
// ------------------------------------Doxygen-Dokumentation---------------------------------------
///  \file csr_helper_functions.hpp
///  \brief
///  Funktionen um das Benutzen des CSR Format zu erleichtern
// ------------------------------------------------------------------------------------------------
// ################################################################################################

#include<iostream>
using namespace std;

// ========================================DOKUMENTATION============================================
//					      csr_format_cpu
///\brief	liefert CSR Format einer Matrix über den Host berechnet
// -------------------------------------------------------------------------------------------------
///@param	matrix    dim x dim Matrix gespeichert als 1D Array der Größe dim*dim
///@param	value	1D Array der Größe "nicht-Null-Einträge", wird bei Laufzeit gefüllt
///@param	column    1D Array der selben Größe wie "value", wird bei Laufzeit gefüllt
///@param pointer   1D Array der Größe "dim+1", wird bei Laufzeit gefüllt
///@param n         Dimension der Matrix
// =================================================================================================

void csr_format_cpu(float* in_matrix,float* value,int* column,int* pointer,int n)
     {
          pointer[0]=0;
          int idx=0;
          for(int i=0;i<n;i++)
          {
               
               for(int j=0;j<n;j++)
               {
                    if(in_matrix[j+i*n]!=0)
                    {
                         value[idx]=in_matrix[j+i*n];
                         column[idx]=j;
                         idx++;
                    }
               }
               pointer[i+1]=idx;     
          }
     }
     
// ========================================DOKUMENTATION============================================
//					      csr_format_cpu
///\brief	liefert CSR Format einer Matrix über den Host berechnet
// -------------------------------------------------------------------------------------------------
///@param	matrix    dim x dim Matrix gespeichert als 1D Array der Größe dim*dim
///@param	n	     Dimension der Matrix

// =================================================================================================    
int nz_values(float* in_matrix,int n)
{
          int count=0;
          for(int i=0;i<n*n;i++)
          {
               if(in_matrix[i]!=0){count++;}
          }
          return count;          
}