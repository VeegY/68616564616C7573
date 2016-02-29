/*
* Projekt:                  Studienprojekt TM 2015/16
* Dateiname:                fullvector.tpp
* Erstellt:                 27.11.15
* Autor / Ansprechpartner:  David
*
* Kurzbeschreibung:
* - Implementiert die Template-Klasse FullVectorGpu und die dazugehörigen
*   VectorTraits.
*/

#ifndef __FULLVECTORGPU_TPP_
#define __FULLVECTORGPU_TPP_

// nur für intellisense
#include "fullvectorgpu.hpp"

namespace Icarus
{

/**********  FullVector  **********/


template<typename Scalar>
FullVectorGpu<Scalar>::FullVectorGpu(size_t dim, MPI_Comm my_comm) :
    _my_comm(my_comm),
    _data(nullptr),
    _dim(dim)
{
    MPI_SCALL(MPI_Comm_rank(_my_comm, &_my_rank));
    MPI_SCALL(MPI_Comm_size(_my_comm, &_num_nodes));
    
    try
    {
        alloc_unifiedV(* _data, _dim);
    }
    catch(...)
    {
        LOG_ERROR("Memory allocation for FullVectorGpu failed.");
    }
}

template<typename Scalar>
FullVectorGpu<Scalar>::FullVectorGpu(const SlicedVector<Scalar>& vec) :
	_my_comm(vec.get_comm()),
    _data(nullptr),
    _dim(vec.get_dim_global())
{
    try
    {
                alloc_unifiedV(* _data, _dim);
    }
	catch (...)
	{
		LOG_ERROR("Memory allocation for FullVectorGpu failed.");
	}

	MPI_SCALL(MPI_Comm_rank(_my_comm, &_my_rank));
	MPI_SCALL(MPI_Comm_size(_my_comm, &_num_nodes));

    Scalar* this_chunk = _data + _my_rank * vec.get_dim_local_nopad();

    // eigenen teil herauskopieren
    for(size_t i=0; i<vec.get_dim_local(); i++)
       this_chunk[i] = vec.get_local(i);
    
	// synchronisiere die teile
    for(int node = 0; node < _num_nodes; node++)
    {
        this_chunk = _data + node * vec.get_dim_local_nopad();
		size_t this_len = (node == _num_nodes - 1) ? vec.get_dim_local_last() : vec.get_dim_local_nopad();
        MPI_SCALL(MPI_Bcast(this_chunk,this_len,ScalarTraits<Scalar>::mpi_type,
                  node,_my_comm));
    }

    // nach der barriere können die kopierquellen gefahrlos überschrieben werden
    // (z.b. matvec-multiplikation mit x = dst)
    MPI_SCALL(MPI_Barrier(_my_comm));
}


template<typename Scalar>
FullVectorGpu<Scalar>::~FullVectorGpu()
{
    if(_data) cleanupgpu(_data);
}

template<typename Scalar>
FullVectorGpu<Scalar>::FullVectorGpu(const FullVectorGpu& other) :
	_my_comm(other._my_comm),
	_my_rank(other._my_rank),
	_num_nodes(other._num_nodes),
    _dim(other._dim)
{
    _data = new Scalar[_dim];
    for (size_t i = 0; i < _dim; i++)
        _data[i] = other._data[i];
}

template<typename Scalar>
FullVectorGpu<Scalar>::FullVectorGpu(FullVectorGpu&& other) :
	_my_comm(other._my_comm),
	_my_rank(other._my_rank),
	_num_nodes(other._num_nodes),
    _dim(other._dim)
{
    _data = other._data;
    other._data = nullptr;
}

template<typename Scalar>
FullVectorGpu<Scalar>&
FullVectorGpu<Scalar>::operator=(const FullVectorGpu& other)
{
    // selbst
    if (this == &other) return *this;
    // fremd
    delete[] _data;
	_my_comm = other._my_comm;
	_my_rank = other._my_rank;
	_num_nodes = other._num_nodes;
    _dim = other._dim;
    
	_data = new Scalar[_dim];
    for (size_t i = 0; i < _dim; i++)
        _data[i] = other._data[i];
    return *this;
}

template<typename Scalar>
FullVectorGpu<Scalar>&
FullVectorGpu<Scalar>::operator=(FullVectorGpu&& other)
{
    // selbst
    if (this == &other) return *this;
    // fremd
	_my_comm = other._my_comm;
	_my_rank = other._my_rank;
	_num_nodes = other._num_nodes;
	_dim = other._dim;

	_data = other._data;
    other._data = nullptr;
    return *this;
}


template<typename Scalar>
typename FullVectorGpu<Scalar>::RealType
FullVectorGpu<Scalar>::l2norm2_impl() const
{
    RealType res(0);
    for(size_t i=0; i<_dim; i++)
        res += ScalarTraits<Scalar>::abs2(_data[i]);
    return res;
}

template<typename Scalar>
typename FullVectorGpu<Scalar>::RealType
FullVectorGpu<Scalar>::maxnorm_impl() const
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
Scalar FullVectorGpu<Scalar>::
scal_prod_impl(const FullVectorGpu<Scalar>& other) const
{
    assert(_dim == other._dim);

    Scalar res(0);
    for(size_t i=0; i<_dim; i++)
        res += ScalarTraits<Scalar>::smult(_data[i], other._data[i]);
    return res;
}

template<typename Scalar>
void FullVectorGpu<Scalar>::
axpy_impl(const Scalar& alpha, const FullVectorGpu<Scalar>& y)
{
    assert(_dim == y._dim);

    for(size_t i=0; i<_dim; i++)
        _data[i] = _data[i] + alpha*y._data[i];
}

template<typename Scalar>
void FullVectorGpu<Scalar>::
scal_impl(const Scalar& alpha)
{
   for(size_t i=0; i<_dim; i++)
        _data[i] *= alpha;
}

template<typename Scalar>
void FullVectorGpu<Scalar>::
copy_impl(const FullVectorGpu<Scalar>& other)
{
   assert(_dim == other._dim);

   for(size_t i=0; i<_dim; i++)
       _data[i] = other._data[i];
}

template<typename Scalar>
void FullVectorGpu<Scalar>::
swap_impl(FullVectorGpu<Scalar>& other)
{
   assert(_dim == other._dim);

   std::swap(_data,other._data);
}


/********** Spezialisierung der VectorTraits **********/


template<typename Scalar>
struct VectorTraits<FullVectorGpu<Scalar>>
{
    typedef typename ScalarTraits<Scalar>::RealType RealType;
    typedef Scalar ScalarType;
};

}

#endif // __FULLVECTORGPU_HPP_
