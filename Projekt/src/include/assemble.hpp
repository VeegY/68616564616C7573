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
	 * \param Nx		Anzahl der (aequidistanten Punkte in x-Richtung)
	 * \param Ny		Anzahl der (aequidistanten Punkte in y-Richtung)
	 * \param Nz		Anzahl der (aequidistanten Punkte in z-Richtung)
	 * \param vertex	Vertex in Node-Zaehlung
	 * \param rhs_val	Wert der rechten Seite am Vertex
	 */
	template<typename Scalar>
	void assemble_row(
		DistEllpackMatrix<Scalar>& A,
		SlicedVector<Scalar>& rhs,
		typename ScalarTraits<Scalar>::RealType h,
		std::vector<char>& types, int Nx, int Ny, size_t vtx_global, Scalar rhs_val);


	/**
	 * \brief	Assembliert eine Matrix und eine rechte Seite
	 *			mittels zentraler Differenzen.
	 *			Zur Zeit werden homogene Dirichlet BCs verwendet
	 *
	 * \param dic_points	Vektor mit dem Typ der auftretenden Punkte
	 * \param h				Schrittweite des finiten Differnezenquotienten
	 * \param Nx			Anzahl der (aequidistanten Punkte in x-Richtung)
	 * \param Ny			Anzahl der (aequidistanten Punkte in y-Richtung)
	 * \param Nz			Anzahl der (aequidistanten Punkte in z-Richtung)
	 * \param rhs_func		Funktion, die jedem raeumlichen Index einen Wert
	 *						der rechten Seite zuweist.
	 *
	 * \return Gibt ein Paar aus Matrix und rechter Seite zurück.
	 */
	template<typename Scalar>
	std::pair < DistEllpackMatrix<Scalar>,
		SlicedVector < Scalar >>
		assemble(std::vector<char>& disc_points,
		typename ScalarTraits<Scalar>::RealType h,
		int Nx, int Ny, int Nz,
		std::function<Scalar(size_t, size_t, size_t)> rhs_func);


	/**
	* \brief	Assembliert eine Matrix	mittels Differenzenquotienten zweiter Ordnung.
	*			Es werden Neumann BC verwendet
	*
	* \param nx			Anzahl der (aequidistanten Punkte in x-Richtung)
	* \param ny			Anzahl der (aequidistanten Punkte in y-Richtung)
	* \param nz			Anzahl der (aequidistanten Punkte in z-Richtung)
	* \param h				Schrittweite des finiten Differnezenquotienten
	* \function bdry			Funktion, die den Neumann-Wert eines Punktes zurückgibt.
	*
	* \return Gibt ein Paar aus Matrix und rechter Seite zurück.
	*/
	template<typename Scalar>
	std::pair <DistEllpackMatrix<Scalar>,
		SlicedVector < Scalar >>
		assemble_neumann(size_t nx, size_t ny, size_t nz,
		typename ScalarTraits<Scalar>::RealType h,
		std::function<Scalar(size_t)> bdry);

    /**
	* \brief	Assembliert eine Matrix	mittels Differenzenquotienten zweiter Ordnung.
	*			Es werden Neumann BC verwendet. Es wird effizienter Assembliert.
	*
	* \param nx			Anzahl der (aequidistanten Punkte in x-Richtung)
	* \param ny			Anzahl der (aequidistanten Punkte in y-Richtung)
	* \param nz			Anzahl der (aequidistanten Punkte in z-Richtung)
	* \param h				Schrittweite des finiten Differnezenquotienten
	* \function bdry			Funktion, die den Neumann-Wert eines Punktes zurückgibt.
	*
	* \return Gibt ein Paar aus Matrix und rechter Seite zurück.
	*/
	template<typename Scalar>
	std::pair <DistEllpackMatrix<Scalar>,
		SlicedVector < Scalar >>
		assemble_neumann_unrolled(size_t nx, size_t ny, size_t nz,
		typename ScalarTraits<Scalar>::RealType h,
		std::function<Scalar(size_t)> bdry);

}

#include "assemble.tpp"

#endif // __ASSEMBLE_HPP_
