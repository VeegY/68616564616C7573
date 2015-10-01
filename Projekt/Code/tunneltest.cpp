#include <iostream>
#include <iomanip>
#include <cmath>
#include "Vector.hpp"
#include "DIA.hpp"
#include "assemble.hpp"
#include <math.h>
#include "PCG.hpp"
#include "vtkFileWriter.hpp"

#define _USE_MATH_DEFINES




int main() 
{	
	 
	
	
	//Gebiet festlegen
	
	int laengex = 3;
	int laengey = 6;
	int laengez = 3;
	
	//gitterweite
	double h = pow(2,-3);
	
	//Gitterpunkte in x - richtung usw
	int Nx = (int) (laengex/h)+1;
	int Ny = (int) (laengey/h)+1;
	int Nz = (int) (laengez/h)+1;
	
	
	
	//algorithmus werte
	Vector<double> T(Nx*Ny*Nz);
	Vector<double> q(T);
	Vector<double> rhs(T);
	Vector<double> c(T);
	Vector<double> k(T);
	int dimD=0;
	
	
	
	
	
	dimD=0;
	//anzahl der Hindernisspunkte bestimmen
	for(int z=0;z<Nz;z++){
		for(int y=0;y<Ny;y++){
			for(int x=0;x<Nx;x++){
				if(x>=Nx/3&&x<=2*Nx/3&&y>=Ny/6&&y<=2*Ny/6&&z>=Nz/3&&z<=2*Nz/3){
					dimD++;
				}
			}
		}
	}
	
	Vector<int> blockid(dimD);
	dimD=0;
	for(int z=0;z<Nz;z++){
		for(int y=0;y<Ny;y++){
			for(int x=0;x<Nx;x++){
				if(x>=Nx/3&&x<=2*Nx/3&&y>=Ny/6&&y<=2*Ny/6&&z>=Nz/3&&z<=2*Nz/3){
					blockid[dimD]=x+Nx*y+Nx*Ny*z;
					dimD++;
				}
			}
		}
	}
	Vector<int> block(Nx*Ny*Nz);
	
	/*Hindernis flag
	1	-	Solid
	0	-	Fluid
	*/
	for(int i=0;i<blockid._dim;i++){
		block[blockid[i]]=1;
	}
	
	
	
	
	//startverteilung setzen
	//Hindernis mit startwaerme
	for(int i=0;i<blockid._dim;i++){
		T[blockid[i]]=100;
	}
	
	

	//konstante fuer luft
	//k.set(0.0262);
	//c.set(1.005);
	
	k.set(0.0262);
	c.set(1.2041*1.005);
	//konstante fuer eisen
	for(int i=0;i<blockid._dim;i++){
		c[blockid[i]]=7874*0.452;
		k[blockid[i]]=80.2;
		q[blockid[i]]=100;
	}
	
	//geschwindigkeiten
	Vector<double> u(rhs);
	Vector<double> v(rhs);
	Vector<double> w(rhs);
	//konstante Geschwindigkeit in y-Richtung
	v.set(100);
	
	//Geschwindigkeit im Hindernis 0 setzen
	for(int z=0;z<Nz;z++){
		for(int y=0;y<Ny;y++){
			for(int x=0;x<Nx;x++){
				if(x>=Nx/3&&x<=2*Nx/3&&y>=Ny/6&&y<=2*Ny/6&&z>=Nz/3&&z<=2*Nz/3){
					v[x+Nx*y+Nx*Ny*z]=0;
					v[x+Nx*(y-1)+Nx*Ny*z]=0;
				}
			}
		}
	}
	
	

	
	//zeitparameter
	//cfl bed: tau<0.5*h*h*roh*c_p/k
	double tau = 0.25*h*h*1.2041*1.005/80.2;
	int Time = 1000;
	vtkFileWriter file ("vtk/tunnel", "test", Nx, Ny, Nz, Time-1);
	
	
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
	int iter = 0;
	
	//start vtk
	file.addPointDataToTimestep(T,iter, "Temperature");
	file.addPointDataToTimestep(block,iter, "Block");
	file.addPointDataToTimestep(v,iter, "y-Velocity");
	iter++;
	
	cout<<"//////////  Begin Algorithm /////////"<<endl;
	for(int time=1;time<Time-1;time++){
		
		//assemble
		assembleT(A, rhs, T, u, v, w, c, k, q, h, tau, Nx, Ny, Nz);
	
		
		//solve
		cout<<PCG_Jacobi(T, A, rhs)<<endl;
		

		//write vtk
		file.addPointDataToTimestep(T,iter, "Temperature");
		file.addPointDataToTimestep(block,iter, "Block");
		file.addPointDataToTimestep(v,iter, "y-Velocity");
		iter++;
		
	}
	cout<<iter<<endl;
	cout<<"///////////   END   ////////////"<<endl;
	return 0;
	

	
}
