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
FullVector<Scalar>::FullVector(size_t dim, MPI_Comm my_comm) :
    _my_comm(my_comm),
    _data(nullptr),
    _dim(dim)
{
    MPI_SCALL(MPI_Comm_rank(_my_comm, &_my_rank));
    MPI_SCALL(MPI_Comm_size(_my_comm, &_num_nodes));
    
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
FullVector<Scalar>::FullVector(const SlicedVector<Scalar>& vec) :
	_my_comm(vec.get_comm()),
    _data(nullptr),
    _dim(vec.get_dim_global())
{
    try
    {
		_data = new Scalar[_dim];
    }
	catch (...)
	{
		LOG_ERROR("Memory allocation for FullVector failed.");
	}

	MPI_SCALL(MPI_Comm_rank(_my_comm, &_my_rank));
	MPI_SCALL(MPI_Comm_size(_my_comm, &_num_nodes));


    Scalar* this_chunk = _data + _my_rank * vec.get_dim_local_nopad();
    
    // eigenen teil herauskopieren
    for(size_t i=0; i<vec.get_dim_local(); i++)
       this_chunk[i] = vec.get_local(i);
    
    // synchronisiere die teile
    if((vec.get_dim_global() % vec.get_dim_local()) != 0)
    {
        LOG_DEBUG("Using non-optimal Bcast version of FullVector construction.");	
        for(int node = 0; node < _num_nodes; node++)
        {
            this_chunk = _data + node * vec.get_dim_local_nopad();
	    	size_t this_len = (node == _num_nodes - 1) ? vec.get_dim_local_last() : vec.get_dim_local_nopad();
            MPI_SCALL(MPI_Bcast(this_chunk,this_len,ScalarTraits<Scalar>::mpi_type,
                      node,_my_comm));
        }
    }
    // Allgather ist schneller, funktioniert aber nur auf gleichen Bloecken.
    else
    {
        MPI_SCALL(MPI_Allgather(this_chunk, vec.get_dim_local(), ScalarTraits<Scalar>::mpi_type, 
          _data, vec.get_dim_local(), ScalarTraits<Scalar>::mpi_type, _my_comm));
    }

    // nach der barriere können die kopierquellen gefahrlos überschrieben werden
    // (z.b. matvec-multiplikation mit x = dst)
    MPI_SCALL(MPI_Barrier(_my_comm));
}


template<typename Scalar>
FullVector<Scalar>::~FullVector()
{
    if(_data) delete[] _data;
}

template<typename Scalar>
FullVector<Scalar>::FullVector(const FullVector& other) :
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
FullVector<Scalar>::FullVector(FullVector&& other) :
	_my_comm(other._my_comm),
	_my_rank(other._my_rank),
	_num_nodes(other._num_nodes),
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
FullVector<Scalar>&
FullVector<Scalar>::operator=(FullVector&& other)
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
