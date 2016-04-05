/*
* Projekt:                  Studienprojekt TM 2015/16
* Dateiname:                fullvector.hpp
* Erstellt:                 27.11.15
* Autor / Ansprechpartner:  David
*
* Kurzbeschreibung:
* - Deklariert die Template-Klasse FullVector.
*/

#ifndef __FULLVECTOR_HPP_
#define __FULLVECTOR_HPP_

#include "logger.hpp"
#include "scalartraits.hpp"
#include <cassert>
#include <memory>
#include <random>
#include <limits>
#include "slicedvector.hpp"
#include "mpihandler.hpp"

namespace Icarus
{

 /**
  * \brief  Vektor, dessen Inhalt komplett auf jeder Node liegt.
  *
  * Alle Elemente dieses Vektors liegen auf jeder Node der zugeordneten Prozessgruppe.
  *
  * \tparam   Scalar  Skalarer Typ der Einträge.
  */
template<typename Scalar>
class FullVector : public Vector<FullVector<Scalar>>
{
    friend class Vector<FullVector<Scalar>>;

    MPI_Comm _my_comm;
    int _my_rank, _num_nodes;

    Scalar* _data;

    size_t _dim;

public:
    typedef Scalar ScalarType;
    typedef typename ScalarTraits<Scalar>::RealType RealType;

   /**
     * \brief   Standardkonstruktor.
     *
     * Erzeugt einen Vektor der Dimension dim, der komplett auf jeder Node der
     * Prozessgruppe my_comm liegt.
     *
     * \param   dim     Dimension des Vektors.
     * \param   my_comm Kommunikator in die Prozessgruppe des Vektors.
     */
    explicit FullVector(size_t dim, MPI_Comm my_comm = MPI_COMM_WORLD);

     /**
     * \brief   Konvertierkonstruktor für einen SlicedVector.
     *
     * Erzeugt einen Vektor mit der globalen Dimension von vec, der alle
     * Teile von vec lokal enthält. Alle Prozesse der Gruppe, die vec verwalten,
     * erhalten eine vollständige Kopie des FullVectors.
     *
     * \param   vec     SlicedVector, der vollständig verteilt werden soll.
     */
    explicit FullVector(const SlicedVector<Scalar>& vec);

    ~FullVector();

    FullVector(const FullVector& other);

    FullVector(FullVector&& other);

    FullVector& operator=(const FullVector& other);

    FullVector& operator=(FullVector&& other);

    /**
     * \brief   Operator für den elementweisen Zugriff.
     *
     * Dieser Operator ermöglicht die elementweise Manipulation des FullVector,
     * analog zu C-Arrays und STL-Containern.
     *
     * \param   index   Index des Elements, auf das zugegriffen werden soll.
     */
    Scalar& operator[](size_t index) {assert(index<_dim); return _data[index];}

    /**
     * \brief   Operator für den elementweisen Zugriff, konstante Variante.
     *
     * Dieser Operator ermöglicht das elementweise Auslesen des FullVector,
     * analog zu C-Arrays und STL-Containern.
     *
     * \param   index   Index des Elements, das gelesen werden soll.
     */
    const Scalar& operator[](size_t index) const { assert(index<_dim); return _data[index]; }

private:
    RealType l2norm2_impl() const;

    RealType maxnorm_impl() const;

    void clear_impl()
    {
        for(size_t i=0; i<_dim; i++) _data[i] = Scalar(0);
    }

    void fill_const_impl(const Scalar& s)
    {
        for(size_t i=0; i<_dim; i++) _data[i] = s;
    }

    Scalar scal_prod_impl(const FullVector& other) const;

    void axpy_impl(const Scalar& alpha, const FullVector& y);

    void scal_impl(const Scalar& alpha);

    void swap_impl(FullVector& other);

    void copy_impl(const FullVector& other);

    size_t get_dim_impl() const { return _dim; }
};

}//namespace Icarus

#include "fullvector.tpp"


#endif // __FULLVECTOR_HPP_
