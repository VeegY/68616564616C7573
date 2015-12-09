#include "distellpackmatrix.hpp"
#include "matrix.hpp"


// 
////////////////////////////////////////////////////////////////////////////////
/*
Assemble Header
Routinen um gegeben Daten entsprechenden Systemmatrizen mit zugehoerigen rechten
Seiten zu assemblieren
Routinen zur Setzung von Randbedingungen
*/
////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////
//Assemblierungsroutinen fuer Poisson Problem (laplace u = 0)
//////////////////////////////////////////////////////////

/*
assmble
////////////////////////////////////
Systemmatrix fuer Poisson Problem 
////////////////////////////////////
assembliert die Zeile der Systemmatrix für einen Knoten und den Eintrag
in der rechte Seite 
aeussere Randzellen werden a priori mit Dirichlet Nullranddaten besetzt
////////////////////////////////////////////////////////////////////////////////
Eingabeparameter
A       -       Systemmatrix	
rhs     -       Rechteseite
h       -       Gitterweite
Type    -       Buchstabe für board/luft/i board
Nx      -       Anzahl Gitterpunkte in x
Ny      -       Anzahl Gitterpunkte in y
Nz      -       Anzahl Gitterpunkte in z Richtung
node    -       Der Aktuelle Knoten für den die Matrix Zeile assembliert werden soll
nodex_  -       linker Nachbarknoten von node in x Richtung
nodex   -       rechter Nachbarknoten von node in x Richtung
nodey_  -       linker Nachbarknoten von node in y Richtung
nodey   -       rechter Nachbarknoten von node in y Richtung
nodez_  -       linker Nachbarknoten von node in z Richtung
nodez   -       rechter Nachbarknoten von node in z Richtung

////////////////////////////////////////////////////////////////////////////////
diskretisiert durch    -     zentrale DQ 2. Ordnung


Knoten benennung startet bei 0
*/

template<typename datatype>
void assemble(datatype* A, datatype* rhs, datatype h,char Type, int Nx, int Ny, int Nz, int node, int nodex_, int nodex, int nodey_, int nodey, int nodez_, int nodez)
 {
   
   const int N = Nx*Ny*Nz; // Dimension der Zeile für node ist 1xN
 
	        // Setze Randwerte
                if (Type=='b'||Type=='o') // Abfrage welcher Typ der Knoten hat
		{
		   
                    A[node]=1.0;

                    rhs[node] = 0.0*(h*h); // Dirichlet Nullranddaten
                }
                else   // Setze innere Punkte durch zentralen DQ 2. Ordnung
                {
		    
		    nodex = node+1;
		    
		    nodey = node+Nx;
		    
		    nodez = node+(Nx*Ny);
		    
		    
                    A[node -(node-nodez)] = 1.0;
                    A[node +(node-nodez)] = 1.0;
                      
                    A[node -(node-nodey)] = 1.0;
                    A[node +(node-nodey)] = 1.0;
                      
                    A[node -(node-nodex)] = 1.0;
                    A[node] = -6.0;
                    A[node +(node-nodex)] = 1.0;
		  
		    rhs[node] = 0.0; // 0 da rechte Seite (f) = 0 ist
                }
                
 
 }