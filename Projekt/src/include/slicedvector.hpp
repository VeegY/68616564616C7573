/*
* Projekt:                  Studienprojekt TM 2015/16
* Dateiname:                slicedvector.hpp
* Erstellt:                 24.11.15
* Autor / Ansprechpartner:  David
*
* Kurzbeschreibung:
* - Deklariert die Template-Klasse SlicedVector.
*/

#ifndef __SLICEDVECTOR_HPP_
#define __SLICEDVECTOR_HPP_

#include "mpihandler.hpp"
#include "logger.hpp"
#include "scalartraits.hpp"
#include <cassert>
#include <memory>
#include <vector>
#include <random>
#include <limits>
#include <ostream>

#include "mpihandler.hpp"
#include "vector.hpp"

namespace Icarus
{

/**
  * \brief  Vektor, dessen Inhalt gleichverteilt auf einer Menge von Nodes liegt.
  *
  * Die Elemente dieses Vektors liegen (annähernd) gleichverteilt auf einer Menge
  * von Nodes der zugeordneten Prozessgruppe.
  *
  * \tparam   Scalar  Skalarer Typ der Einträge.
  */
template<typename Scalar>
class SlicedVector : public Vector<SlicedVector<Scalar>>
{
    friend class Vector<SlicedVector<Scalar>>;

	// mpi umgebung
	MPI_Comm _my_comm;
	int _my_rank, _num_nodes;

    // Globale und lokale Dimension
    size_t _dim_global, _dim_local;
    
    // Größe des lokalen Blocks auf den ersten N-1 bzw. der letzten Node
    size_t _dim_local_nopad, _dim_local_last;
    Scalar* _data;

public:
    typedef Scalar ScalarType;
    typedef typename ScalarTraits<Scalar>::RealType RealType;

   /**
     * \brief   Standardkonstruktor.
     *
     * Konstuiert einen Vektor, dessen Elemente auf den Nodes der Prozessgruppe
     * my_comm verteilt werden. Die ersten N-1 Nodes enthalten die gleiche Anzahl
     * an Elementen, die letzte Node enthält den (eventuell etwas kleineren) Rest.
     *
     * \param   dim_global   Globale Dimension des Vektors, also Summe der Größen der
     *                       auf die Nodes verteilten Blöcke.
     */
    explicit SlicedVector(size_t dim_global, MPI_Comm my_comm = MPI_COMM_WORLD);

    ~SlicedVector();

    SlicedVector(const SlicedVector& other);

    SlicedVector(SlicedVector&& other);

    SlicedVector& operator=(const SlicedVector& other);

    SlicedVector& operator=(SlicedVector&& other);

    /**
      * \brief   Setze den Wert val an die globale Position pos.
      *
      * Setze den Wert val an die globale Position pos. Diese Operation erfordert 
      * MPI-Kommunikation, wenn die zu setzende Position nicht auf der Node liegt,
      * die den Befehl ausführt, und ist daher hinsichtlich Effizienz mit Vorsicht
      * zu benutzen.
      *
      * \param   pos   Globale Position des Elements, das gesetzt werden soll.
      * \param   val   Wert, der an die Stelle pos kopiert werden soll.
      */    
    void set_global(size_t pos, const Scalar& val);

    /**
      * \brief   Hole den Eintrag an der globalen Position pos.
      *
      * Hole den Wert an der globale Position pos. Diese Operation erfordert 
      * MPI-Kommunikation, wenn die zu lesende Position nicht auf der Node liegt,
      * die den Befehl ausführt, und ist daher hinsichtlich Effizienz mit Vorsicht
      * zu benutzen.
      *
      * \param   pos   Globale Position des Elements, das gelesen werden soll.
      *
      * \return Eintrag an der globalen Position pos.
      */
    Scalar get_global(size_t pos) const;

    /**
      * \brief   Gibt die globale Dimension des Vektors zurück.
      */
    size_t get_dim_global() const {return _dim_global;}

    /**
      * \brief   Gibt die lokale Dimension des Vektors, d.h. die Größe
      *          des auf der aufrufenden Node gespeicherten Blocks zurück.
      */
    size_t get_dim_local() const {return _dim_local;}

    /**
      * \brief   Gibt die Größe eines Blocks, wie er auf den ersten N-1 Nodes liegt, zurück.
      */
    size_t get_dim_local_nopad() const {return _dim_local_nopad;}

    /**
      * \brief   Gibt die Größe des Blocks, der auf der letzten Node liegt, zurück.
      */
	size_t get_dim_local_last() const { return _dim_local_last; }

    /**
      * \brief   Gibt den MPI-Kommunikator in doe Prozessgrupe, der der Vektor gehört, zurück.
      */
	MPI_Comm get_comm() const { return _my_comm; }

    /**
      * \brief   Schreibe den lokalen Inhalt des Block in den Stream out.
      *
      * Für die Verwendung dieser Funktion muss eine entsprechende Überladung des
      * Operators std::ostream::operator<<(Scalar) existieren.
      *
      * \param out  Stream, in den die Ausgabe geschrieben werden soll.
      */
	void print_local_data(std::ostream& out) const;

    /**
      * \brief   Setze den Wert val an die lokale Position pos.
      *
      * Setze den Wert val an die lokale Position pos. Diese Operation erfordert 
      * keine MPI-Kommunikation.
      *
      * \param   pos   Lokale Position des Elements, das gesetzt werden soll.
      * \param   val   Wert, der an die Stelle pos kopiert werden soll.
      */ 
    void set_local(size_t pos, const Scalar& val)
    {
        _data[pos] = val;
    }

    /**
      * \brief   Hole den Eintrag an der lokalen Position pos.
      *
      * Hole den Wert an der lokalen Position pos. Diese Operation erfordert 
      * keine MPI-Kommunikation.
      *
      * \param   pos   Lokale Position des Elements, das gelesen werden soll.
      *
      * \return Eintrag an der lokalen Position pos.
      */
    Scalar get_local(size_t pos) const
    {
        return _data[pos];
    }

private:
    void swap_impl(SlicedVector& other);

    void copy_impl(const SlicedVector& other);

    RealType l2norm2_impl() const;

    RealType maxnorm_impl() const;

    void clear_impl()
    {
        for(size_t i=0; i<_dim_local; i++) _data[i] = Scalar(0);
    }

    void fill_const_impl(const Scalar& s)
    {
        for(size_t i=0; i<_dim_local; i++) _data[i] = s;
    }

    Scalar scal_prod_impl(const SlicedVector& other) const;

    void axpy_impl(const Scalar& alpha, const SlicedVector& y);

    void scal_impl(const Scalar& alpha);

    size_t get_dim_impl() const { return _dim_global; }

};

}

#include "slicedvector.tpp"

#endif // __SLICEDVECTOR_H_
