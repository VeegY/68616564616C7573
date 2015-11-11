#include <iostream>
#include <iomanip>
#include <cmath>
#include "vector.hpp"
#include "dia.hpp"
#include "assemble.hpp"
#include <math.h>
#include "pcg.hpp"
#include "vtkfilewriter.hpp"

#define _USE_MATH_DEFINES




int main() 
{	
	 
	
	
	//Gebiet festlegen
	
	int laengex = 1;
	int laengey = 1;
	int laengez = 1;
	
	//gitterweite
	double h = pow(2,-4);
	
	//Gitterpunkte in x - richtung usw
	//hier alle richtungen gleich
	int Nx = (int) (laengex/h)+1;
	int Ny = Nx;
	int Nz = Nx;
	
	
	
	//algorithmus werte
	Vector<double> T(Nx*Ny*Nz);
	Vector<double> q(T);
	Vector<double> rhs(T);
	
	
	//startverteilung setzen
	for(int z=0;z<Nz;z++){
		for(int y=0;y<Ny;y++){
			for(int x=0;x<Nx;x++){
				T[x+ Nx*y+ Nx*Ny*z] = sin(M_PI*x*h)*sin(M_PI*y*h)*sin(M_PI*z*h);
			}
		}
	}
	
	//loesung
	Vector<double> sol(T);
	Vector<double> error(T);
	int k=0;
	
	//zeitparameter
	double tau = 0.25*h*h;
	int Time = 100;
	
	vtkFileWriter file ("vtk/heattest1", "test", Nx, Ny, Nz, Time-1);
	
	//Matrix initialisierung
	Vector<double> daten(7*Nx*Ny*Nz);
	  Vector<int> offset(7);
	  offset[0]=-Nx*Ny;
	  offset[1]=-Nx;
	  offset[2]=-1;
	  offset[3]=-0;
	  offset[4]=1;
	  offset[5]=Nx;
	  offset[6]=Nx*Ny;
	DIA<double> A(Nx*Ny*Nz,7,daten,offset);
	
	
	
	//werte fuer Loeser
	double TOL = 0.0001;
	double norm = 0;
	int iter = 0;
	Vector<double> x(T);
	for(int j=0;j<T._dim;j++){
		error[j]=sqrt(pow((T[j]-sol[j]*exp(-3*M_PI*M_PI*tau*iter)),2));
			
	}
	
	//write start-vtk
	file.addPointDataToTimestep(T,iter, "Temperature");
	file.addPointDataToTimestep(error,iter, "error");
	iter++;
	
	cout<<"//////////  Begin Algorithm /////////"<<endl;
	for(int time=1;time<Time-1;time++){
	
	
		assembleT(A,rhs,T,q,h,tau,Nx,Ny,Nz);
		
		//solve
		PCG_Jacobi(T, A, rhs);
		
		
		for(int j=0;j<T._dim;j++){
			norm+=pow((T[j]-sol[j]*exp(-3*M_PI*M_PI*tau*iter)),2);
			error[j]=sqrt(pow((T[j]-sol[j]*exp(-3*M_PI*M_PI*tau*iter)),2));
			
		}
		
		if(sqrt(norm)/Nx/Ny/Nz>0.001){	//h=2^-4;
			cout<<"Test - failed"<<endl;
			return -1;
		}
		
		norm=0;
		//write vtk
		file.addPointDataToTimestep(T,iter, "Temperature");
		file.addPointDataToTimestep(error,iter, "error");
		iter++;
		
	}
	cout<<iter<<endl;
	
	cout<<"///////////   END   ////////////"<<endl;
	return 0;
	

	
}
