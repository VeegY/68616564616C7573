#ifndef __ASSEMBLE_HPP_
#define __ASSEMBLE_HPP_

#include "distellpackmatrix.hpp"
#include "slicedvector.hpp"
#include "scalartraits.hpp"

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

namespace Icarus
{

template<typename Scalar, int _num_nodes, int _first_node = 0>
void assemble(
        DistEllpackMatrix<Scalar,_num_nodes,_first_node>& A,
        SlicedVector<Scalar, _num_nodes, _first_node>& rhs,
        ScalarTraits<Scalar>::RealType h,
        char Type, int Nx, int Ny, int Nz, int node)
{
    A.prepare_sequential_fill(7);

    // Setze Randwerte
    if (Type=='b'||Type=='o') // Abfrage welcher Typ der Knoten hat
    {
        A.sequential_fill(node, 1.0);
        rhs.set_local(node,0.0*(h*h)); // Dirichlet Nullranddaten
    }

    else   // Setze innere Punkte durch zentralen DQ 2. Ordnung
    {
        int nodex = node+1;
        int nodey = node+Nx;
        int nodez = node+(Nx*Ny);

        A.sequential_fill(node -(node-nodez), 1.0);
        A.sequential_fill(node +(node-nodez), 1.0);

        A.sequential_fill(node -(node-nodey), 1.0);
        A.sequential_fill(node +(node-nodey), 1.0);

        A.sequential_fill(node -(node-nodex), 1.0);
        A.sequential_fill(node, -6.0);
        A.sequential_fill(node +(node-nodex),1.0);

        rhs.set_local(node,0.0); // 0 da rechte Seite (f) = 0 ist
    }
}

}
#endif // __ASSEMBLE_HPP_
