#include <mpi.h>
#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <iomanip>
using namespace std;

void Fibonacci(int* fib);


int main(int argc, char* argv[])
{

  MPI_Init(&argc, &argv);
 
  int rank;
  int tag = 99;
  
   
  MPI_Status status;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
//RANK1//====================================================
  
  if(rank == 0)
  {  
     int *fib = new int [2];
     fib[0]=0;
     fib[1]=1;  
      
     MPI_Send(fib,
             2,
             MPI_INT,
             1,
             tag,
             MPI_COMM_WORLD);
     
     
     MPI_Recv(fib,
             2,
             MPI_INT,
             11,
             tag,
             MPI_COMM_WORLD,
             &status);
             
             cout << fib[0] << " " << fib[1] << "\n";
  }
  
//RANK2//====================================================   
     else if(rank >= 1 && rank <=10)
     {
          int targetbuffer[2];
      
          MPI_Recv(targetbuffer,
             2,
             MPI_INT,
             rank-1,
             tag,
             MPI_COMM_WORLD,
             &status);
             
          Fibonacci(targetbuffer);
          
          MPI_Send(targetbuffer,
             2,
             MPI_INT,
             rank+1,
             tag,
             MPI_COMM_WORLD);                             
     }
     
     else if(rank == 11)
     {
          int targetbuffer[2];
          
          MPI_Recv(targetbuffer,
             2,
             MPI_INT,
             rank-1,
             tag,
             MPI_COMM_WORLD,
             &status);
             
          Fibonacci(targetbuffer);
          
          MPI_Send(targetbuffer,
             2,
             MPI_INT,
             0,
             tag,
             MPI_COMM_WORLD); 
     
     }

  MPI_Finalize(); 
  return 0;
}
