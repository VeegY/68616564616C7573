#ifndef __ASSEMBLE_HPP_
#define __ASSEMBLE_HPP_

#include "distellpackmatrix.hpp"
#include "slicedvector.hpp"
#include "scalartraits.hpp"
#include "utility.hpp"

namespace Icarus
{

/**
 * \brief	Assembliert eine Zeile (ein Vertex) mittels zentraler Differenzen.
 *			Zur Zeit werden homogene Dirichlet BCs verwendet
 *
 * \param A			Matrix, die aufgebaut wird
 * \param rhs		Rechte Seite, die aufgebaut wird
 * \param h			Schrittweite des finiten Differnezenquotienten
 * \param Type		Typ des Vertex
 * \param Nx		Anzahl der (äquidistanten Punkte in x-Richtung)
 * \param Ny		Anzahl der (äquidistanten Punkte in y-Richtung)
 * \param Nz		Anzahl der (äquidistanten Punkte in z-Richtung)
 * \param vertex	Vertex in Node-Zählung
 * \param rhs_val	Wert der rechten Seite am Vertex
 */
template<typename Scalar, int _num_nodes, int _first_node = 0>
void assemble_row(
        DistEllpackMatrix<Scalar,_num_nodes,_first_node>& A,
        SlicedVector<Scalar, _num_nodes, _first_node>& rhs,
        typename ScalarTraits<Scalar>::RealType h,
        char Type, int Nx, int Ny, int vertex, Scalar rhs_val)
{
    // Setze Randwerte
    if (Type=='b'||Type=='o') // Abfrage welcher Typ der Knoten hat
    {
        A.sequential_fill(vertex, 1.0);
        rhs.set_local(vertex,0.0*(h*h)); // Dirichlet Nullranddaten
    }

    else   // Setze innere Punkte durch zentralen DQ 2. Ordnung
    {
        int vertexx = vertex+1;
        int vertexy = vertex+Nx;
        int vertexz = vertex+(Nx*Ny);

        A.sequential_fill(vertex -(vertex-vertexz), 1.0);
        A.sequential_fill(vertex +(vertex-vertexz), 1.0);

        A.sequential_fill(vertex -(vertex-vertexy), 1.0);
        A.sequential_fill(vertex +(vertex-vertexy), 1.0);

        A.sequential_fill(vertex -(vertex-vertexx), 1.0);
        A.sequential_fill(vertex, -6.0);
        A.sequential_fill(vertex +(vertex-vertexx),1.0);

        rhs.set_local(vertex,rhs_val); // 0 da rechte Seite (f) = 0 ist
    }
}

/**
 * \brief	Assembliert eine Matrix und eine rechte Seite
 *			mittels zentraler Differenzen.
 *			Zur Zeit werden homogene Dirichlet BCs verwendet
 *
 * \param dic_points	Vektor mit dem Typ der auftretenden Punkte
 * \param h				Schrittweite des finiten Differnezenquotienten
 * \param Nx			Anzahl der (äquidistanten Punkte in x-Richtung)
 * \param Ny			Anzahl der (äquidistanten Punkte in y-Richtung)
 * \param Nz			Anzahl der (äquidistanten Punkte in z-Richtung)
 * \param rhs_func		Funktion, die jedem räumlichen Index einen Wert
 *						der rechten Seite zuweist.
 *
 * \return Gibt ein Paar aus Matrix und rechter Seite zurück.
 */
template<typename Scalar, int _num_nodes, int _first_node = 0>
std::pair<DistEllpackMatrix<Scalar,_num_nodes,_first_node>,
SlicedVector<Scalar, _num_nodes, _first_node>>
assemble(std::vector<char>& disc_points,
        typename ScalarTraits<Scalar>::RealType h,
        int Nx, int Ny, int Nz,
		std::function<Scalar(size_t,size_t,size_t)> rhs_func)
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

	A.prepare_sequential_fill(7);

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
                assemble_row(A,rhs,h,disc_points[index],Nx,Ny,index-fron,rhs_func(x,y,z));
            }
        }
    }
    return {A, rhs};
}

}
#endif // __ASSEMBLE_HPP_
