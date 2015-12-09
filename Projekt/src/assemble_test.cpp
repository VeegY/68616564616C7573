#include "include/assemble2.hpp"
#include <iostream> // Für Test ausgaben
#include<sstream>
#include<fstream>

int main()
{
///// Dimensionen Übergeben //////


const  long int Nx(10),Ny(Nx),Nz(Nx); // Setze Dimensionen in x,y,z - Richtung
const long int N =Nx*Ny*Nz;

double h = 1.0/(Nx+1); // Schrittweite h


// Abfrage ob genug speicher Angelegt werden kann
  try
{
  double *A = new double[N]; // Matrix der Dimension Nx*Ny*N
  double *rhs = new double[N]; // rechte Seite  dim N
}
catch(...)
{
  std::cout<<"Allocation failed "<< std::endl;
  return 0;
}
// Lege Speicher an
double *A = new double[N]; // Matrix der Dimension Nx*Ny*N
double *rhs = new double[N]; // rechte Seite  dim N



 char Eigen ='a';
   assemble(A,rhs,h,Eigen,Nx,Ny,Nz,200,0,0,0,0,0,0); // Assembliere Zeile für den i-ten Knoten




// Ausgabe zum testen //
/*
 std::fstream fout;
 fout.open("Matrix.txt");


 for (int i(0);i<N;++i)
 {

   assemble2(A,rhs,h,Eigen,Nx,Ny,Nz,i,0,0,0,0,0,0); // Assembliere Zeile für den i-ten Knoten
   for (int spal(0); spal<N; spal++)
	      {

		fout << A[spal];
	      }
   for (int j(0);j<N;++j)
   {
     A[j] = 0;
   }
   fout << std::endl;
 }

 fout.close();

 */
delete[] A,rhs;
}
