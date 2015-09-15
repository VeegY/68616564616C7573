#include <iostream>
#include <iomanip>
#include <cmath>
#include "Vector.hpp"
#include "DIA.hpp"
#include "assemble.hpp"
#include <math.h>

#define _USE_MATH_DEFINES




int main() 
{	
	//Gebiet festlegen
	
	int laengex = 1;
	int laengey = 1;
	int laengez = 1;
	
	//gitterweite
	double h = pow(2,-3);
	
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
	
	//zeitparameter
	double tau = 0.1;
	int Time = 15;
	
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
	
	//write start-vtk
	/*...*/
	iter++;
	
	cout<<"//////////  Begin Algorithm /////////"<<endl;
	for(int time=0;time<Time;time++){
	
	
		assembleT(A,rhs,T,q,h,tau,Nx,Ny,Nz);
		//solve
		for(int m =0;m<100; m++){
			for(int i=0;i<T._dim;i++){
				for(int j=0;j<T._dim;j++){
					if(i!=j){
						x[i]=(rhs[i]-A.value(i,j)*T[j])/A.value(i,i);
					}
				}
			}
			T=x;
		}
		
		for(int j=0;j<T._dim;j++){
			norm+=(T[j]-sol[j]*exp(-M_PI*M_PI*tau*iter))*(T[j]-sol[j]*exp(-M_PI*M_PI*tau*iter));
		}
		//cout<<norm;
		norm=0;
		
		//write vtk
		/*...*/
		iter++;
		
	}
	
	
	cout<<"///////////   END   ////////////"<<endl;
	
	

	
}
