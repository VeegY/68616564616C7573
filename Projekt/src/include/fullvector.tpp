/*
* Projekt:                  Studienprojekt TM 2015/16
* Dateiname:                fullvector.tpp
* Erstellt:                 27.11.15
* Autor / Ansprechpartner:  David
*
* Kurzbeschreibung:
* - Implementiert die Template-Klasse FullVector und die dazugehörigen
*   VectorTraits.
*/

#ifndef __FULLVECTOR_TPP_
#define __FULLVECTOR_TPP_

// nur für intellisense
#include "fullvector.hpp"

namespace Icarus
{

/**********  FullVector  **********/


template<typename Scalar>
FullVector<Scalar>::FullVector(size_t dim) :
    _data(nullptr),
    _dim(dim)
{
    try
    {
        _data = new Scalar[_dim];
    }
    catch(...)
    {
        LOG_ERROR("Memory allocation for FullVector failed.");
    }
}

template<typename Scalar>
template<int _num_nodes, int _first_node>
FullVector<Scalar>::FullVector(const SlicedVector<Scalar, _num_nodes, _first_node>& vec) :
    _dim(vec.get_dim_global()),
    _data(nullptr)
{
    try
    {
        _dim = vec.get_dim_global();
        _data = new Scalar[_dim];
    }
    catch(...)
    {
        LOG_ERROR("Memory allocation for FullVector failed.");
    }

    int this_node = MPI_HANDLER.get_my_rank();
    Scalar* this_chunk = _data + (this_node - _first_node) * vec.get_dim_local_nopad();

    // haben wir selbst einen teil des slicedvectors?
    if(this_node >= _first_node && this_node < _first_node + _num_nodes)
    {
        // dann herauskopieren
        for(size_t i=0; i<vec.get_dim_local(); i++)
            this_chunk[i] = vec.get_local(i);
    }

    // synchronisiere die teile
    for(int node = _first_node; node < _first_node + _num_nodes; node++)
    {
        this_chunk = _data + (node - _first_node) * vec.get_dim_local_nopad();
        MPI_SCALL(MPI_Bcast(this_chunk,vec.get_dim_local(),ScalarTraits<Scalar>::mpi_type,
                  node,MPI_COMM_WORLD));
    }

    // nach der barriere können die kopierquellen gefahrlos überschrieben werden
    // (z.b. matvec-multiplikation mit x = dst)
    MPI_SCALL(MPI_Barrier(MPI_COMM_WORLD));
}


template<typename Scalar>
FullVector<Scalar>::~FullVector()
{
    if(_data) delete[] _data;
}

template<typename Scalar>
FullVector<Scalar>::FullVector(const FullVector& other) :
    _dim(other._dim)
{
    _data = new Scalar[_dim];
    for (size_t i = 0; i < _dim; i++)
        _data[i] = other._data[i];
}

template<typename Scalar>
FullVector<Scalar>::FullVector(const FullVector&& other) :
    _dim(other._dim)
{
    _data = other._data;
    other._data = nullptr;
}

template<typename Scalar>
FullVector<Scalar>&
FullVector<Scalar>::operator=(const FullVector& other)
{
    // selbst
    if (this == &other) return *this;
    // fremd
    delete[] _data;
    _dim = other._dim;
    _data = new Scalar[_dim];
    for (size_t i = 0; i < _dim; i++)
        _data[i] = other._data[i];
    return *this;
}

template<typename Scalar>
FullVector<Scalar>&
FullVector<Scalar>::operator=(const FullVector&& other)
{
    // selbst
    if (this == &other) return *this;
    // fremd
    _dim = other._dim;
    _data = other._data;
    other._data = nullptr;
    return *this;
}


template<typename Scalar>
typename FullVector<Scalar>::RealType
FullVector<Scalar>::l2norm2_impl() const
{
    RealType res(0);
    for(size_t i=0; i<_dim; i++)
        res += ScalarTraits<Scalar>::abs2(_data[i]);
    return res;
}

template<typename Scalar>
typename FullVector<Scalar>::RealType
FullVector<Scalar>::maxnorm_impl() const
{
    RealType res = std::numeric_limits<RealType>::min();
    for(size_t i=0; i<_dim; i++)
    {
        RealType tmp = ScalarTraits<Scalar>::abs(_data[i]);
        if(tmp > res) res = tmp;
    }
    return res;
}

template<typename Scalar>
Scalar FullVector<Scalar>::
scal_prod_impl(const FullVector<Scalar>& other) const
{
    assert(_dim == other._dim);

    Scalar res(0);
    for(size_t i=0; i<_dim; i++)
        res += ScalarTraits<Scalar>::smult(_data[i], other._data[i]);
    return res;
}

template<typename Scalar>
void FullVector<Scalar>::
axpy_impl(const Scalar& alpha, const FullVector<Scalar>& y)
{
    assert(_dim == y._dim);

    for(size_t i=0; i<_dim; i++)
        _data[i] = _data[i] + alpha*y._data[i];
}

template<typename Scalar>
void FullVector<Scalar>::
scal_impl(const Scalar& alpha)
{
   for(size_t i=0; i<_dim; i++)
        _data[i] *= alpha;
}

template<typename Scalar>
void FullVector<Scalar>::
copy_impl(const FullVector<Scalar>& other)
{
   assert(_dim == other._dim);

   for(size_t i=0; i<_dim; i++)
       _data[i] = other._data[i];
}

template<typename Scalar>
void FullVector<Scalar>::
swap_impl(FullVector<Scalar>& other)
{
   assert(_dim == other._dim);

   std::swap(_data,other._data);
}


/********** Spezialisierung der VectorTraits **********/


template<typename Scalar>
struct VectorTraits<FullVector<Scalar>>
{
    typedef typename ScalarTraits<Scalar>::RealType RealType;
    typedef Scalar ScalarType;
};

}

#endif // __FULLVECTOR_HPP_
