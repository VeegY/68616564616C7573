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
#include "vector.hpp"

namespace Icarus
{

template<typename Scalar, int _num_nodes, int _first_node = 0>
class SlicedVector : public Vector<SlicedVector<Scalar, _num_nodes, _first_node>>
{
    friend class Vector<SlicedVector<Scalar, _num_nodes, _first_node>>;

    static const int _last_node = _first_node + _num_nodes - 1;

    size_t _dim_global, _dim_local, _dim_local_nopad;
    Scalar* _data;

public:
    typedef Scalar ScalarType;
    typedef typename ScalarTraits<Scalar>::RealType RealType;

    explicit SlicedVector(size_t dim_global);

    ~SlicedVector();

    explicit SlicedVector(const SlicedVector& other);

    explicit SlicedVector(const SlicedVector&& other);

    SlicedVector& operator=(const SlicedVector& other);

    SlicedVector& operator=(const SlicedVector&& other);

    void set_global(size_t pos, const Scalar& val);

    Scalar get_global(size_t pos) const;

    size_t get_dim_global() const {return _dim_global;}

    size_t get_dim_local() const {return _dim_local;}

    size_t get_dim_local_nopad() const {return _dim_local_nopad;}

    void set_local(size_t pos, const Scalar& val)
    {
        _data[pos] = val;
    }

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
