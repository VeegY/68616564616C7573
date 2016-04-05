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


template<typename Scalar>
void gpu_ax_(Scalar *data, Scalar *fvec, Scalar *result, int *indices, int max_row_length, int dim_local);

template<typename type>
void gpu_dot_(type *vecx, type *vecy, size_t dim, type *erg);

template<typename type>
void gpu_axpy(type *vecx, type scalar, type *vecy, size_t dim);

template<typename type>
void gpu_l2(type *vec, size_t dim, type *erg);

template<typename type>
void gpumaxnorm(type *vec, size_t dim, type *erg);

template<typename type>
void copygpu_(type *vecin, type *vecout, size_t dim);

template <typename Scalar>
void cleanupgpu(Scalar *data);

template<typename Scalar>
void alloc_unified(Scalar **fvec, size_t dim_fvec);

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
        alloc_unified(& _data, _dim);
    }
    catch(...)
    {
        LOG_ERROR("Memory allocation for FullVectorGpu failed.");
    }
}

template<typename Scalar>
FullVectorGpu<Scalar>::FullVectorGpu(const SlicedVectorGpu<Scalar>& vec) :
	_my_comm(vec.get_comm()),
    _data(nullptr),
    _dim(vec.get_dim_global())
{
    try
    {
        alloc_unified(& _data, _dim);
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
    try
    {
        alloc_unified(& _data, _dim);
    }
    catch(...)
    {
        LOG_ERROR("Memory allocation for SlicedVectorGpu failed.");
    }
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
    if (_data) cleanupgpu(_data);
	_my_comm = other._my_comm;
	_my_rank = other._my_rank;
	_num_nodes = other._num_nodes;
    _dim = other._dim;

	try
    {
        alloc_unified(& _data, _dim);
    }
    catch(...)
    {
        LOG_ERROR("Memory allocation for SlicedVectorGpu failed.");
    }
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
    RealType *res, res2;
    alloc_unified(&res, 1);

    gpu_l2(_data,_dim, res);
    res2=*res;
    cleanupgpu(res);
    return res2;
}


template<typename Scalar>
typename FullVectorGpu<Scalar>::RealType
FullVectorGpu<Scalar>::maxnorm_impl() const
{
    RealType *res, res2;
    alloc_unified(&res, 1);
    *res = std::numeric_limits<RealType>::min();

    gpumaxnorm(_data,_dim, res);
    res2 = *res;
    cleanupgpu(res);
    return res2;
}

template<typename Scalar>
Scalar FullVectorGpu<Scalar>::
scal_prod_impl(const FullVectorGpu<Scalar>& other) const
{
    assert(_dim == other._dim);

    Scalar *erg, erg2;
    alloc_unified(&erg, 1.0);

    gpu_dot_(_data, other.getDataPointer(), _dim, erg);
    erg2 = *erg;
    cleanupgpu(erg);
    return erg2;
}


template<typename Scalar>
void FullVectorGpu<Scalar>::
axpy_impl(const Scalar& alpha, const FullVectorGpu<Scalar>& y)
{
    assert(_dim == y._dim);

    //Scalar alpha2(alpha);
    FullVectorGpu<Scalar> yvec(y);
    //alloc_unified((Scalar **)&alpha2, (size_t)1.0);

    gpu_axpy(yvec.getDataPointer(), alpha, _data, _dim);

   // cleanupgpu(&alpha2);
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

   copygpu_(other.getDataPointer(), _data,  _dim);
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
