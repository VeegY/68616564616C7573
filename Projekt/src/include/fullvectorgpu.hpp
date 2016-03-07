/*
* Projekt:                  Studienprojekt TM 2015/16
* Dateiname:                fullvector.hpp
* Erstellt:                 27.11.15
* Autor / Ansprechpartner:  David
*
* Kurzbeschreibung:
* - Deklariert die Template-Klasse FullVector.
*/

#ifndef __FULLVECTORGPU_HPP_
#define __FULLVECTORGPU_HPP_

#include "logger.hpp"
#include "scalartraits.hpp"
#include <cassert>
#include <memory>
#include <random>
#include <limits>
#include "slicedvectorgpu.hpp"
#include "mpihandler.hpp"
#include "distellpackmatrixgpu.hpp"

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
class FullVectorGpu : public Vector<FullVectorGpu<Scalar>>
{
    friend class Vector<FullVectorGpu<Scalar>>;

    MPI_Comm _my_comm;
    int _my_rank, _num_nodes;

    Scalar* _data;

    size_t _dim;

public:
    typedef Scalar ScalarType;
    typedef typename ScalarTraits<Scalar>::RealType RealType;

  //  friend DistEllpackMatrixGpu<Scalar>::mult_vec_impl(const FullVectorGpu&, SlicedVectorGpu&);
    
   /**
     * \brief   Standardkonstruktor.
     *
     * Erzeugt einen Vektor der Dimension dim, der komplett auf jeder Node der
     * Prozessgruppe my_comm liegt.
     *
     * \param   dim     Dimension des Vektors.
     * \param   my_comm Kommunikator in die Prozessgruppe des Vektors.
     */
    explicit FullVectorGpu(size_t dim, MPI_Comm my_comm = MPI_COMM_WORLD);

     /**
     * \brief   Konvertierkonstruktor für einen SlicedVectorGpu.
     *
     * Erzeugt einen Vektor mit der globalen Dimension von vec, der alle
     * Teile von vec lokal enthält. Alle Prozesse der Gruppe, die vec verwalten,
     * erhalten eine vollständige Kopie des FullVectorGpus.
     *
     * \param   vec     SlicedVectorGpu, der vollständig verteilt werden soll.
     */
    explicit FullVectorGpu(const SlicedVectorGpu<Scalar>& vec);

    ~FullVectorGpu();

    FullVectorGpu(const FullVectorGpu& other);

    FullVectorGpu(FullVectorGpu&& other);

    FullVectorGpu& operator=(const FullVectorGpu& other);

    FullVectorGpu& operator=(FullVectorGpu&& other);
    
    /**
     * \brief   Operator für den elementweisen Zugriff.
     *
     * Dieser Operator ermöglicht die elementweise Manipulation des FullVectorGpu,
     * analog zu C-Arrays und STL-Containern.
     *
     * \param   index   Index des Elements, auf das zugegriffen werden soll.
     */
    Scalar& operator[](size_t index) {assert(index<_dim); return _data[index];}

    /**
     * \brief   Operator für den elementweisen Zugriff, konstante Variante.
     *
     * Dieser Operator ermöglicht das elementweise Auslesen des FullVectorGpu,
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

    Scalar scal_prod_impl(const FullVectorGpu& other) const;

    void axpy_impl(const Scalar& alpha, const FullVectorGpu& y);

    void scal_impl(const Scalar& alpha);

    void swap_impl(FullVectorGpu& other);

    void copy_impl(const FullVectorGpu& other);

    size_t get_dim_impl() const { return _dim; }
};

}

#include "fullvectorgpu.tpp"


#endif // __FULLVECTORGPU_HPP_
