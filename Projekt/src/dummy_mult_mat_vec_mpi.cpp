#include <mpi.h>
#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <iomanip>
using namespace std;

template<typename Scalar>
void mult_vec_memcpy(Scalar *data, Scalar *fvec, Scalar *result, int *indices, int max_row_length, int dim_local);


int main(int argc, char* argv[])
{

  //MPI_Init(&argc, &argv);
 
  //int rank;
  //int tag = 99;
  
   
  //MPI_Status status;
  //MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
//RANK1//====================================================
  
  //if(rank == 0)
  //{  
     int max_row_length = 3;
     int dim_local = 5;
     
     
     int data_i[15]{1,2,3,4,5,0,6,7,0,8,0,0,9,0,0};
     float data_f[15]{1.0,2.0,3.0,4.0,5.0,0.0,6.0,7.0,0.0,8.0,0.0,0.0,9.0,0.0,0.0};
     
     int fvec_i[5]{1,2,3,4,5};
     float fvec_f[5]{1.0,2.0,3.0,4.0,5.0};
     
     int result_i[5];
     float result_f[5];
     
     
     int indices[15]{2,4,5,1,4,0,1,4,0,1,0,0,2,0,0};
     
     
     mult_vec_memcpy<int>(data_i, fvec_i, result_i, indices, max_row_length, dim_local);
     mult_vec_memcpy<float>(data_f, fvec_f, result_f, indices, max_row_length, dim_local);
     
     
     
     
        
    
 // }

    

  //MPI_Finalize(); 
  return 0;
}
