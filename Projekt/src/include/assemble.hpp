#ifndef __ASSEMBLE_HPP_
#define __ASSEMBLE_HPP_

#include "distellpackmatrix.hpp"
#include "slicedvector.hpp"
#include "scalartraits.hpp"
#include "utility.hpp"

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
void assemble_row(
        DistEllpackMatrix<Scalar,_num_nodes,_first_node>& A,
        SlicedVector<Scalar, _num_nodes, _first_node>& rhs,
        typename ScalarTraits<Scalar>::RealType h,
        char Type, int Nx, int Ny, int node)
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

template<typename Scalar, int _num_nodes, int _first_node = 0>
std::pair<DistEllpackMatrix<Scalar,_num_nodes,_first_node>,
SlicedVector<Scalar, _num_nodes, _first_node>>
assemble(std::vector<char>& disc_points,
        typename ScalarTraits<Scalar>::RealType h,
        int Nx, int Ny, int Nz)
{
    const size_t N = Nx*Ny*Nz;
    DistEllpackMatrix<Scalar,_num_nodes,_first_node> A(N);
    SlicedVector<Scalar, _num_nodes, _first_node> rhs(N);

    const size_t fron = A.first_row_on_node();
    size_t fron_x, fron_y, fron_z;
    deflatten_3d(fron, Nx, Ny, fron_x, fron_y, fron_z);

    const size_t end = fron + A.get_dim_local();
    size_t end_x, end_y, end_z;
    deflatten_3d(end, Nx, Ny, end_x, end_y, end_z);

    for(size_t z=fron_z; z<end_z; z++)
    {
        size_t ymin = (z==fron_z) ? fron_y : 0;
        size_t ymax = (z==end_z-1) ? end_y : Ny;
        for(size_t y=ymin; y<ymax; y++)
        {
            size_t xmin = (z==fron_z && y==fron_y) ? fron_x : 0;
            size_t xmax = (z==end_z && y==end_y-1) ? end_x : Nx;
            for(size_t x=xmin; x<xmax; x++)
            {
                const size_t index = x + y*Nx + z*Nx*Ny;
                assemble_row(A,rhs,h,disc_points[index],Nx,Ny,index-fron);
            }
        }
    }
    return {A, rhs};
}

}
#endif // __ASSEMBLE_HPP_
