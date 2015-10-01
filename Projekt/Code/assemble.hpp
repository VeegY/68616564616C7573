#include <iostream>
#include <iomanip>
#include <cmath>
#include "Vector.hpp"
#include "DIA.hpp"
#include <algorithm>


////////////////////////////////////////////////////////////////////////////////
/*
Assemble Header

Routinen um gegeben Daten entsprechenden Systemmatrizen mit zugehoerigen rechten
Seiten zu assemblieren
Zeitschrittverfahren: Impl. Euler

Routinen zur Setzung von Randbedingungen
*/
////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////
//Assemblierungsroutinen fuer Transport/Waermeleitungsproblem
//////////////////////////////////////////////////////////

/*
assmbleT
////////////////////////////////////
Systemmatrix fuer Transportgleichung 
////////////////////////////////////

assembliert neue Systemmatrix und rechte Seite und ueberschreibt alte Werte
Datenfelden muessen bei Eingabe richtige Groesse haben

aeussere Randzellen werden a priori mit Dirichlet Nullranddaten besetzt
////////////////////////////////////////////////////////////////////////////////
Eingabeparameter
A	-	Systemmatrix	
rhs	-	Rechteseite
T	-	Temperatur
u	-	Geschwindigkeit in x
v	-	Geschwindigkeit in y
w	-	Geschwindigkeit in z Richtung
c	-	spezifische Waermekapazitaet*Dichte
k	-	Waermeausbreitungskoeefizient
q	-	Waermequelle
h	-	Gitterweite
tau	-	Zeitschrittweite
Nx	-	Anzahl Gitterpunkte in x
Ny	-	Anzahl Gitterpunkte in y
Nz	-	Anzahl Gitterpunkte in z Richtung
////////////////////////////////////////////////////////////////////////////////

Diffusionsterm	-	zentrale DQ 2. Ordnung
Konvektion	-	Upwind Scheme

*/

template<typename data>
void assembleT(DIA<data>& A, Vector<data>& rhs, Vector<data>& T, Vector<data>& u, Vector<data>& v, Vector<data>& w, Vector<data>& c, Vector<data>& k, Vector<data>& q, double h, double tau, int Nx, int Ny, int Nz)
{	
	
	int N = Nx*Ny*Nz;
	for(int z=0;z<Nz;z++){
		for(int y=0;y<Ny;y++){
			for(int x=0;x<Nx;x++){
				//Randzellen
				if(x==0||y==0||z==0||x==Nx-1||y==Ny-1||z==Nz-1){
					(*(A._data))[x+Nx*y+Nx*Ny*z]=0;
					(*(A._data))[x+Nx*y+Nx*Ny*z+N]=0;
					(*(A._data))[x+Nx*y+Nx*Ny*z+2*N]=0;
					(*(A._data))[x+Nx*y+Nx*Ny*z+3*N]=1;
					(*(A._data))[x+Nx*y+Nx*Ny*z+4*N]=0;
					(*(A._data))[x+Nx*y+Nx*Ny*z+5*N]=0;
					(*(A._data))[x+Nx*y+Nx*Ny*z+6*N]=0;
					//Rechte seite setzen
					rhs[x+Nx*y+Nx*Ny*z]=0;
				}
				//innere Punkte
				else{
					//setzen der werte beginnend bei der untersten nebendiagonalen
					(*(A._data))[x+Nx*y+Nx*Ny*z]=-c[x+Nx*y+Nx*Ny*z]/k[x+Nx*y+Nx*Ny*z]/h/h\
									-max(w[x+Nx*y+Nx*Ny*z],0.0)/h;
					(*(A._data))[x+Nx*y+Nx*Ny*z+N]=-c[x+Nx*y+Nx*Ny*z]/k[x+Nx*y+Nx*Ny*z]/h/h\
									-max(v[x+Nx*y+Nx*Ny*z],0.0)/h;
					(*(A._data))[x+Nx*y+Nx*Ny*z+2*N]=-c[x+Nx*y+Nx*Ny*z]/k[x+Nx*y+Nx*Ny*z]/h/h
									-max(u[x+Nx*y+Nx*Ny*z],0.0)/h;;
					(*(A._data))[x+Nx*y+Nx*Ny*z+3*N]=6*c[x+Nx*y+Nx*Ny*z]/k[x+Nx*y+Nx*Ny*z]/h/h\
									+max(u[x+Nx*y+Nx*Ny*z],0.0)/h\
									+max(v[x+Nx*y+Nx*Ny*z],0.0)/h\
									+max(w[x+Nx*y+Nx*Ny*z],0.0)/h\
									-min(u[x+Nx*y+Nx*Ny*z],0.0)/h\
									-min(v[x+Nx*y+Nx*Ny*z],0.0)/h\
									-min(w[x+Nx*y+Nx*Ny*z],0.0)/h\
									+1/tau;
					(*(A._data))[x+Nx*y+Nx*Ny*z+4*N]=-c[x+Nx*y+Nx*Ny*z]/k[x+Nx*y+Nx*Ny*z]/h/h\
									+min(u[x+Nx*y+Nx*Ny*z],0.0)/h;
					(*(A._data))[x+Nx*y+Nx*Ny*z+5*N]=-c[x+Nx*y+Nx*Ny*z]/k[x+Nx*y+Nx*Ny*z]/h/h\
									+min(v[x+Nx*y+Nx*Ny*z],0.0)/h;
					(*(A._data))[x+Nx*y+Nx*Ny*z+6*N]=-c[x+Nx*y+Nx*Ny*z]/k[x+Nx*y+Nx*Ny*z]/h/h\
									+min(w[x+Nx*y+Nx*Ny*z],0.0)/h;
					//Rechte seite setzen
					rhs[x+Nx*y+Nx*Ny*z]=T[x+Nx*y+Nx*Ny*z]/tau+q[x+Nx*y+Nx*Ny*z]/c[x+Nx*y+Nx*Ny*z];
				}
			}
		}
	}
	
}

/*
assmbleT
////////////////////////////
Waermeleitungsgleichung mit waermeuebertragungskoeffizietn c ungleich const 1
////////////////////////////
assembliert neue Systemmatrix und rechte Seite und ueberschreibt alte Werte
Datenfelden muessen bei Eingabe richtige Groesse haben

aeussere Randzellen werden a priori mit Dirichlet Nullranddaten besetzt
////////////////////////////////////////////////////////////////////////////////
Eingabeparameter
A
rhs
T	-	Temperatur
c	-	Waermeuebertragungskonstante
q	-	Waermequelle
h	-	Gitterweite
tau	-	Zeitschrittweite
Nx	-	Anzahl Gitterpunkte in x
Ny	-	Anzahl Gitterpunkte in y
Nz	-	Anzahl Gitterpunkte in z Richtung
////////////////////////////////////////////////////////////////////////////////

Diffusionsterm	-	zentrale DQ 2. Ordnung


*/
template<typename data>
void assembleT(DIA<data>& A, Vector<data>& rhs, Vector<data>& T, Vector<data>& c, Vector<data>& q, double h, double tau, int Nx, int Ny, int Nz)
{	
	
	int N = Nx*Ny*Nz;
	for(int z=0;z<Nz;z++){
		for(int y=0;y<Ny;y++){
			for(int x=0;x<Nx;x++){
				//Randzellen
				if(x==0||y==0||z==0||x==Nx-1||y==Ny-1||z==Nz-1){
					(*(A._data))[x+Nx*y+Nx*Ny*z]=0;
					(*(A._data))[x+Nx*y+Nx*Ny*z+N]=0;
					(*(A._data))[x+Nx*y+Nx*Ny*z+2*N]=0;
					(*(A._data))[x+Nx*y+Nx*Ny*z+3*N]=1;
					(*(A._data))[x+Nx*y+Nx*Ny*z+4*N]=0;
					(*(A._data))[x+Nx*y+Nx*Ny*z+5*N]=0;
					(*(A._data))[x+Nx*y+Nx*Ny*z+6*N]=0;
					//Rechte seite setzen
					rhs[x+Nx*y+Nx*Ny*z]=0;
				}
				//innere Punkte
				else{
					//setzen der werte beginnend bei der untersten nebendiagonalen
					(*(A._data))[x+Nx*y+Nx*Ny*z]=-c[x+Nx*y+Nx*Ny*z]/h/h;
					(*(A._data))[x+Nx*y+Nx*Ny*z+N]=-c[x+Nx*y+Nx*Ny*z]/h/h;
					(*(A._data))[x+Nx*y+Nx*Ny*z+2*N]=-c[x+Nx*y+Nx*Ny*z]/h/h;
					(*(A._data))[x+Nx*y+Nx*Ny*z+3*N]=6*c[x+Nx*y+Nx*Ny*z]/h/h\
									+1/tau;
					(*(A._data))[x+Nx*y+Nx*Ny*z+4*N]=-c[x+Nx*y+Nx*Ny*z]/h/h;
					(*(A._data))[x+Nx*y+Nx*Ny*z+5*N]=-c[x+Nx*y+Nx*Ny*z]/h/h;
					(*(A._data))[x+Nx*y+Nx*Ny*z+6*N]=-c[x+Nx*y+Nx*Ny*z]/h/h;
					//Rechte seite setzen
					rhs[x+Nx*y+Nx*Ny*z]=T[x+Nx*y+Nx*Ny*z]/tau+q[x+Nx*y+Nx*Ny*z];
				}
			}
		}
	}
	
}



/*
assmbleT
////////////////////////////
Waermeleitungsgleichung
////////////////////////////
assembliert neue Systemmatrix und rechte Seite und ueberschreibt alte Werte
Datenfelden muessen bei Eingabe richtige Groesse haben

aeussere Randzellen werden a priori mit Dirichlet Nullranddaten besetzt
////////////////////////////////////////////////////////////////////////////////
Eingabeparameter
A
rhs
T	-	Temperatur
q	-	Waermequelle
h	-	Gitterweite
tau	-	Zeitschrittweite
Nx	-	Anzahl Gitterpunkte in x
Ny	-	Anzahl Gitterpunkte in y
Nz	-	Anzahl Gitterpunkte in z Richtung
////////////////////////////////////////////////////////////////////////////////

Diffuisionsterm	-	zentrale DQ 2. Ordnung	

*/
template<typename data>
void assembleT(DIA<data>& A, Vector<data>& rhs, Vector<data>& T, Vector<data>& q, double h, double tau, int Nx, int Ny, int Nz)
{	
	
	
	int N = Nx*Ny*Nz;
	for(int z=0;z<Nz;z++){
		for(int y=0;y<Ny;y++){
			for(int x=0;x<Nx;x++){
				//Randzellen
				if(x==0||y==0||z==0||x==Nx-1||y==Ny-1||z==Nz-1){
					(*(A._data))[x+Nx*y+Nx*Ny*z]=0;
					(*(A._data))[x+Nx*y+Nx*Ny*z+N]=0;
					(*(A._data))[x+Nx*y+Nx*Ny*z+2*N]=0;
					(*(A._data))[x+Nx*y+Nx*Ny*z+3*N]=1;
					(*(A._data))[x+Nx*y+Nx*Ny*z+4*N]=0;
					(*(A._data))[x+Nx*y+Nx*Ny*z+5*N]=0;
					(*(A._data))[x+Nx*y+Nx*Ny*z+6*N]=0;
					//Rechte seite setzen
					rhs[x+Nx*y+Nx*Ny*z]=0;
				}
				//innere Punkte
				else{
					//setzen der werte beginnend bei der untersten nebendiagonalen
					(*(A._data))[x+Nx*y+Nx*Ny*z]=-1/h/h;
					(*(A._data))[x+Nx*y+Nx*Ny*z+N]=-1/h/h;
					(*(A._data))[x+Nx*y+Nx*Ny*z+2*N]=-1/h/h;
					(*(A._data))[x+Nx*y+Nx*Ny*z+3*N]=6/h/h\
									+1/tau;
					(*(A._data))[x+Nx*y+Nx*Ny*z+4*N]=-1/h/h;
					(*(A._data))[x+Nx*y+Nx*Ny*z+5*N]=-1/h/h;
					(*(A._data))[x+Nx*y+Nx*Ny*z+6*N]=-1/h/h;
					//Rechte seite setzen
					rhs[x+Nx*y+Nx*Ny*z]=T[x+Nx*y+Nx*Ny*z]/tau+q[x+Nx*y+Nx*Ny*z];
				}
				
			}
		}
	}
	
}


//////////////////////////////////////////////////
/*
introduceContraints
Randsetzungsroutinen - Problemunabhängig
*/
//////////////////////////////////////////////////



////////////
/*
TO DO:
Neumannrandaten in der Matrix setzen
*/
////////////

/*
introduceConstraints
///////////////////////
mit Dirichlet und Neumannrandwerten
///////////////////////
modifiziert Systemmatrix und rechte Seite fuer Randwerte, dabei wird Matrixzeile
entsprechend veraendert und der entsprechende Wert in die rechte Seite geschrieben
////////////////////////////////////////////////////////////////////////////////
Eingabeparameter:
A	-	Systemmatrix
rhs	-	rechte Seite
h	-	Gitterweite
indexD	-	Gitterindizes in denen Dirichlet randdaten vorgegben werden
dirichlet-	Dirichlet Randdaten
indexN	-	Gitterindizes in denen Neumann randdaten vorgegben werden
normals	-	Normalenrichtung
neumann	-	Neumann Randdaten

////////////////////////////////////////////////////////////////////////////////
*/
template<typename data>
void introduceConstraints(DIA<data>& A, Vector<data>& rhs, double h, Vector<int>& indexD, Vector<data>& dirichlet,Vector<int>& indexN, Vector<int>& normals, Vector<data>& neumann)
{	
	
	//Dirichlet Randdaten
	for(int d=0;d<indexD._dim;d++){
		//einheitszeilen setzen
		for(int i=0;i<A._numDiags;i++){
			if((*(A._offset))[i]==0){
				(*(A._data))[indexD[d]+i*A._dim]=1;
			}
			else{
				(*(A._data))[indexD[d]+i*A._dim]=0;
			}
		}
		//rhs modifizieren
		rhs[indexD[d]]=dirichlet[d];
	}
	//Neumann Randdaten
	for(int n=0;n<indexN._dim;n++){
		for(int i=0;i<A._numDiags;i++){
			//abfragen ueber normalen Richtung
			/*...*/
		}
		//rhs modifizieren
		rhs[indexN[n]]=neumann[n];
	}
}


/*
introduceConstraints
////////////////////////
Nur Dirichlet Randwerte
////////////////////////
modifiziert Systemmatrix und rechte Seite fuer Randwerte, dabei wird Matrixzeile
entsprechend veraendert und der entsprechende Wert in die rechte Seite geschrieben
////////////////////////////////////////////////////////////////////////////////
Eingabeparameter:
A	-	Systemmatrix
rhs	-	rechte Seite
h	-	Gitterweite
indexD	-	Gitterindizes in denen Dirichlet randdaten vorgegben werden
dirichlet-	Dirichlet Randdaten

////////////////////////////////////////////////////////////////////////////////
*/
template<typename data>
void introduceConstraints(DIA<data>& A, Vector<data>& rhs, double h, Vector<int>& indexD, Vector<data>& dirichlet)
{	

	//Dirichlet Randdaten
	for(int d=0;d<indexD._dim;d++){
		//einheitszeilen setzen
		for(int i=0;i<A._numDiags;i++){
			if((*(A._offset))[i]==0){
				(*(A._data))[indexD[d]+i*A._dim]=1;
			}
			else{
				(*(A._data))[indexD[d]+i*A._dim]=0;
			}
		}
		//rhs modifizieren
		rhs[indexD[d]]=dirichlet[d];
	}
}


////////////
/*
TO DO:
Neumannrandaten in der Matrix setzen
*/
////////////
/*
introduceConstraints
///////////////////////
mit Dirichlet und Neumannrandwerten
///////////////////////
modifiziert Systemmatrix und rechte Seite fuer Randwerte, dabei wird Matrixzeile
entsprechend veraendert und der entsprechende Wert in die rechte Seite geschrieben
////////////////////////////////////////////////////////////////////////////////
Eingabeparameter:
A	-	Systemmatrix
rhs	-	rechte Seite
h	-	Gitterweite
indexN	-	Gitterindizes in denen Neumann randdaten vorgegben werden
normals	-	Normalenrichtung
neumann	-	Neumann Randdaten

////////////////////////////////////////////////////////////////////////////////
*/
template<typename data>
void introduceConstraints(DIA<data>& A, Vector<data>& rhs, double h, Vector<int>& indexN, Vector<int>& normals, Vector<data>& neumann)
{	
	//Neumann Randdaten
	for(int n=0;n<indexN._dim;n++){
		for(int i=0;i<A._numDiags;i++){
			//abfragen ueber normalen Richtung
			/*...*/
		}
		//rhs modifizieren
		rhs[indexN[n]]=neumann[n];
	}
}




