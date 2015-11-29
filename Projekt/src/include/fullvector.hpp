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

namespace Icarus
{

template<typename Scalar>
class FullVector : public Vector<FullVector<Scalar>>
{
    friend class Vector<FullVector<Scalar>>;

    Scalar* _data;

    size_t _dim;

public:
    typedef Scalar ScalarType;
    typedef typename ScalarTraits<Scalar>::RealType RealType;

    explicit FullVector(size_t dim);

    template<int _num_nodes, int _first_node>
    explicit FullVector(const SlicedVector<Scalar, _num_nodes, _first_node>& vec);

    ~FullVector();

    explicit FullVector(const FullVector& other);

    explicit FullVector(const FullVector&& other);

    FullVector& operator=(const FullVector& other);

    FullVector& operator=(const FullVector&& other);

    Scalar& operator[](size_t index) {assert(index<_dim); return _data[index];}

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

}

#include "fullvector.tpp"


#endif // __FULLVECTOR_HPP_
