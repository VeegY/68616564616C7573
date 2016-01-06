#ifndef __ASSEMBLE_HPP_
#define __ASSEMBLE_HPP_

#include <functional>

#include "mpihandler.hpp"
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
template<typename Scalar>
void assemble_row(
        DistEllpackMatrix<Scalar>& A,
        SlicedVector<Scalar>& rhs,
        typename ScalarTraits<Scalar>::RealType h,
		std::vector<char>& types, int Nx, int Ny, size_t vtx_global, Scalar rhs_val)
{
	const size_t fron = A.first_row_on_node();
	const size_t vtx_local = vtx_global - fron;

	// Typ boundary or object
    if (types[vtx_global] == 'b' || types[vtx_global] == 'o')
    {
        A.sequential_fill(vtx_local, 1.0);
		// dirichlet
        rhs.set_local(vtx_local,rhs_val);
    }
	// Typ freier knoten
    else  
    {
		// nachbarn (x+,x-,y+,y-,z+,z-)
		const size_t nn[6] = {
			vtx_global + 1, vtx_global - 1,
			vtx_global + Nx, vtx_global - Nx,
			vtx_global + Nx*Ny, vtx_global - Nx*Ny };
        
		A.sequential_fill(vtx_local, -6.0);
		for (int i = 0; i < 6; i++) A.sequential_fill(nn[i], 1.0);

		// rechte seite
        rhs.set_local(vtx_local,rhs_val);
    }
	A.end_of_row();
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
template<typename Scalar>
std::pair<DistEllpackMatrix<Scalar>,
SlicedVector<Scalar>>
assemble(std::vector<char>& disc_points,
        typename ScalarTraits<Scalar>::RealType h,
        int Nx, int Ny, int Nz,
		std::function<Scalar(size_t,size_t,size_t)> rhs_func)
{
    const size_t N = Nx*Ny*Nz;
    DistEllpackMatrix<Scalar> A(N);
    SlicedVector<Scalar> rhs(N);

    const size_t fron = A.first_row_on_node();
    size_t fron_x, fron_y, fron_z;
    deflatten_3d(fron, Nx, Ny, fron_x, fron_y, fron_z);

    const size_t lron = fron + A.get_dim_local() - 1;
    size_t lron_x, lron_y, lron_z;
    deflatten_3d(lron, Nx, Ny, lron_x, lron_y, lron_z);

	// 7 punkte stern
	const unsigned max_nnz_per_line = 7;
	A.prepare_sequential_fill(max_nnz_per_line);

    for(size_t z=fron_z; z<=lron_z; z++)
    {
        size_t ymin = (z==fron_z) ? fron_y : 0;
        size_t ymax = (z==lron_z) ? lron_y : Ny-1;
        for(size_t y=ymin; y<=ymax; y++)
        {
            size_t xmin = (z==fron_z && y==fron_y) ? fron_x : 0;
            size_t xmax = (z==lron_z && y==lron_y) ? lron_x : Nx-1;
            for(size_t x=xmin; x<=xmax; x++)
            {
                const size_t index = x + y*Nx + z*Nx*Ny;
                assemble_row(A,rhs,h,disc_points,Nx,Ny,index,rhs_func(x,y,z));
				//std::cout << "Assembling row " << index << " (" << x << "," << y << "," << z <<") on node " << MPI_HANDLER.get_my_rank() << std::endl;
			}
        }
    }
    return {A, rhs};
}

}
#endif // __ASSEMBLE_HPP_
