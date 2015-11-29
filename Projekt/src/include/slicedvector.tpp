/*
* Projekt:                  Studienprojekt TM 2015/16
* Dateiname:                slicedvector.tpp
* Erstellt:                 24.11.15
* Autor / Ansprechpartner:  David
*
* Kurzbeschreibung:
* - Implementiert die Template-Klasse SlicedVector und die dazugehörigen
*   VectorTraits.
*/

#ifndef __SLICEDVECTOR_TPP_
#define __SLICEDVECTOR_TPP_

// nur für intellisense
#include "slicedvector.hpp"

namespace Icarus
{

/**********  SlicedVector  **********/


template<typename Scalar, int _num_nodes, int _first_node>
SlicedVector<Scalar, _num_nodes, _first_node>::
SlicedVector(size_t dim_global) :
    _dim_global(dim_global),
    _dim_local(0),
    _dim_local_nopad(0),
    _data(nullptr)
{
    // wenn weniger nodes vorhanden als angefordert, abbruch
    if(_last_node + 1 > MPI_HANDLER.get_n_procs())
        LOG_ERROR("SlicedVector is not compatible with node structure.");

    // ist diese node überhaupt beteiligt?
    if(MPI_HANDLER.get_my_rank() > _last_node || MPI_HANDLER.get_my_rank() < _first_node) return;

    // geht die division genau auf
    if(_dim_global % _num_nodes == 0)
    {
        _dim_local = _dim_local_nopad = _dim_global/_num_nodes;
    }
    // funktioniert nur, wenn dim_global >> _num_nodes
    else
    {
        _dim_local = _dim_local_nopad = _dim_global/_num_nodes + 1;
        if(MPI_HANDLER.get_my_rank() == _last_node)
            _dim_local = _dim_global - (_num_nodes - 1)*_dim_local_nopad;
        assert(_dim_local >= 0);
    }
    // allokiere lokalen speicher
    try
    {
        _data = new Scalar[_dim_local];
    }
    catch(...)
    {
        LOG_ERROR("Memory allocation for SlicedVector failed.");
    }
    // LOG_DEBUG("dim_local ",_dim_local,", dim_local_nopad ", _dim_local_nopad);
}

template<typename Scalar, int _num_nodes, int _first_node>
SlicedVector<Scalar, _num_nodes, _first_node>::
~SlicedVector()
{
    if(_data) delete[] _data;
}

template<typename Scalar, int _num_nodes, int _first_node>
SlicedVector<Scalar, _num_nodes, _first_node>::
SlicedVector(const SlicedVector& other) :
    _dim_global(other._dim_global),
    _dim_local(other._dim_local),
    _dim_local_nopad(other._dim_local_nopad)
{
    _data = new Scalar[_dim_local];
    for (size_t i = 0; i < _dim_local; i++)
        _data[i] = other._data[i];
}

template<typename Scalar, int _num_nodes, int _first_node>
SlicedVector<Scalar, _num_nodes, _first_node>::
SlicedVector(const SlicedVector&& other) :
    _dim_global(other._dim_global),
    _dim_local(other._dim_local),
    _dim_local_nopad(other._dim_local_nopad)
{
    _data = other._data;
    other._data = nullptr;
}

template<typename Scalar, int _num_nodes, int _first_node>
SlicedVector<Scalar, _num_nodes, _first_node>&
SlicedVector<Scalar, _num_nodes, _first_node>::
operator=(const SlicedVector& other)
{
    // selbst
    if (this == &other) return this;
    // fremd
    delete[] _data;
    _dim_global = other._dim_global;
    _dim_local = other._dim_local;
    _dim_local_nopad = other._dim_local_nopad;
    _data = new Scalar[_dim_local];
    for (size_t i = 0; i < _dim_local; i++)
        _data[i] = other._data[i];
    return *this;
}

template<typename Scalar, int _num_nodes, int _first_node>
SlicedVector<Scalar, _num_nodes, _first_node>&
SlicedVector<Scalar, _num_nodes, _first_node>::
operator=(const SlicedVector&& other)
{
    // selbst
    if (this == &other) return this;
    // fremd
    _dim_global = other._dim_global;
    _dim_local = other._dim_local;
    _dim_local_nopad = other._dim_local_nopad;
    _data = other._data;
    other._data = nullptr;
    return *this;
}

template<typename Scalar, int _num_nodes, int _first_node>
void SlicedVector<Scalar, _num_nodes, _first_node>::
set_global(size_t pos, const Scalar& val)
{
    int affected_rank = pos / _dim_local_nopad;
    if (MPI_HANDLER.get_my_rank() == affected_rank)
        _data[pos - (affected_rank-_first_node)*_dim_local_nopad] = val;
}

template<typename Scalar, int _num_nodes, int _first_node>
Scalar SlicedVector<Scalar, _num_nodes, _first_node>::
get_global(size_t pos) const
{
    int affected_rank = pos / _dim_local_nopad;
    Scalar val;
    if (MPI_HANDLER.get_my_rank() == affected_rank)
    {
        val = _data[pos - (affected_rank - _first_node)*_dim_local_nopad];
    }
    MPI_SCALL(MPI_Bcast(&val,1,ScalarTraits<Scalar>::mpi_type,
                        affected_rank,MPI_COMM_WORLD));
    return val;
}

template<typename Scalar, int _num_nodes, int _first_node>
typename SlicedVector<Scalar, _num_nodes, _first_node>::RealType
SlicedVector<Scalar, _num_nodes, _first_node>::
l2norm2_impl() const
{
    RealType res(0), res_global;
    for(size_t i=0; i<_dim_local; i++)
        res += ScalarTraits<Scalar>::abs2(_data[i]);
    MPI_SCALL(MPI_Allreduce(&res, &res_global, 1,
                            ScalarTraits<RealType>::mpi_type, MPI_SUM, MPI_COMM_WORLD));
    return res_global;
}

template<typename Scalar, int _num_nodes, int _first_node>
typename SlicedVector<Scalar, _num_nodes, _first_node>::RealType
SlicedVector<Scalar, _num_nodes, _first_node>::
maxnorm_impl() const
{
    RealType res = std::numeric_limits<RealType>::min(), res_global, tmp;
    for(size_t i=0; i<_dim_local; i++)
    {
        tmp = ScalarTraits<Scalar>::abs(_data[i]);
        if(tmp > res) res = tmp;
    }
    MPI_SCALL(MPI_Allreduce(&res, &res_global, 1,
                            ScalarTraits<RealType>::mpi_type, MPI_MAX, MPI_COMM_WORLD));
    return res_global;
}

template<typename Scalar, int _num_nodes, int _first_node>
Scalar SlicedVector<Scalar,_num_nodes, _first_node>::
scal_prod_impl(const SlicedVector<Scalar, _num_nodes, _first_node>& other) const
{
    assert(_dim_global == other._dim_global);

    Scalar res(0), res_global;
    for(size_t i=0; i<_dim_local; i++)
        res += ScalarTraits<Scalar>::smult(_data[i], other._data[i]);
    MPI_SCALL(MPI_Allreduce(&res, &res_global, 1,
                            ScalarTraits<Scalar>::mpi_type, MPI_SUM, MPI_COMM_WORLD));
    return res_global;
}

template<typename Scalar, int _num_nodes, int _first_node>
void SlicedVector<Scalar, _num_nodes, _first_node>::
axpy_impl(const Scalar& alpha, const SlicedVector<Scalar, _num_nodes, _first_node>& y)
{
    assert(_dim_global == y._dim_global);

    for(size_t i=0; i<_dim_local; i++)
        _data[i] = _data[i] + alpha*y._data[i];
}

template<typename Scalar, int _num_nodes, int _first_node>
void SlicedVector<Scalar, _num_nodes, _first_node>::
scal_impl(const Scalar& alpha)
{
   for(size_t i=0; i<_dim_local; i++)
        _data[i] *= alpha;
}


template<typename Scalar, int _num_nodes, int _first_node>
void SlicedVector<Scalar, _num_nodes, _first_node>::
swap_impl(SlicedVector &other)
{
    assert(_dim_global == other._dim_global);

    std::swap(_data,other._data);
}

template<typename Scalar, int _num_nodes, int _first_node>
void SlicedVector<Scalar, _num_nodes, _first_node>::
copy_impl(const SlicedVector &other)
{
    assert(_dim_global == other._dim_global);

    for(size_t i=0; i<_dim_local; i++) _data[i] = other._data[i];
}


/********** Spezialisierung der VectorTraits **********/


template<typename Scalar, int _num_nodes, int _first_node>
struct VectorTraits<SlicedVector<Scalar,_num_nodes,_first_node>>
{
    typedef typename ScalarTraits<Scalar>::RealType RealType;
    typedef Scalar ScalarType;
};

}

#endif // __SLICEDVECTOR_TPP_
