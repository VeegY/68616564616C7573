/*
* Projekt:                  Studienprojekt TM 2015/16
* Dateiname:                slicedvector.tpp
* Erstellt:                 24.11.15
* Autor / Ansprechpartner:  David
*
* Kurzbeschreibung:
* - Implementiert die Template-Klasse slicedvectorgpu und die dazugehörigen
*   VectorTraits.
*/

#ifndef __SLICEDVECTORGPU_TPP_
#define __SLICEDVECTORGPU_TPP_

// nur für intellisense
#include "slicedvectorgpu.hpp"


template<typename Scalar>
void gpu_ax_(Scalar *data, Scalar *fvec, Scalar *result, int *indices, int max_row_length, int dim_local);

template <typename Scalar>
void cleanupgpu(Scalar *data);

template<typename Scalar>
void alloc_unifiedD(Scalar **data, int **indices, int max_row_length, int dim_local);

template<typename Scalar>
void alloc_unifiedV(Scalar **fvec, int dim_fvec);

namespace Icarus
{

/**********  slicedvectorgpu  **********/


template<typename Scalar>
SlicedVectorGpu<Scalar>::
SlicedVectorGpu(size_t dim_global, MPI_Comm my_comm) :
	_my_comm(my_comm),
    _dim_global(dim_global),
    _dim_local(0),
    _dim_local_nopad(0),
	_dim_local_last(0),
    _data(nullptr)
{
	// hole informationen über die mpi umgebung
	MPI_SCALL(MPI_Comm_size(_my_comm, &_num_nodes));
	MPI_SCALL(MPI_Comm_rank(_my_comm, &_my_rank));

	// geht die division genau auf
    if(_dim_global % _num_nodes == 0)
    {
        _dim_local = _dim_local_nopad = _dim_local_last = _dim_global/_num_nodes;
    }
    // funktioniert nur, wenn dim_global >> _num_nodes
    else
    {
        _dim_local = _dim_local_nopad = _dim_global/_num_nodes + 1;
		_dim_local_last = _dim_global - (_num_nodes - 1)*_dim_local_nopad;
        if(_my_rank == _num_nodes - 1)
            _dim_local = _dim_local_last;
        assert(_dim_local >= 0);
    }
    // allokiere lokalen speicher
    try
    {
        alloc_unifiedV(& _data, _dim_local);
    }
    catch(...)
    {
        LOG_ERROR("Memory allocation for SlicedVectorGpu failed.");
    }
    // LOG_DEBUG("dim_local ",_dim_local,", dim_local_nopad ", _dim_local_nopad);
}

template<typename Scalar>
SlicedVectorGpu<Scalar>::
~SlicedVectorGpu()
{
    if(_data) cleanupgpu(_data);
}

template<typename Scalar>
SlicedVectorGpu<Scalar>::
SlicedVectorGpu(const SlicedVectorGpu& other) :
	_my_comm(other._my_comm),
	_my_rank(other._my_rank),
	_num_nodes(other._num_nodes),
    _dim_global(other._dim_global),
    _dim_local(other._dim_local),
    _dim_local_nopad(other._dim_local_nopad),
	_dim_local_last(other._dim_local_last)
{
    _data = new Scalar[_dim_local];
    for (size_t i = 0; i < _dim_local; i++)
        _data[i] = other._data[i];
}

template<typename Scalar>
SlicedVectorGpu<Scalar>::
SlicedVectorGpu(SlicedVectorGpu&& other) :
	_my_comm(other._my_comm),
	_my_rank(other._my_rank),
	_num_nodes(other._num_nodes),
    _dim_global(other._dim_global),
    _dim_local(other._dim_local),
    _dim_local_nopad(other._dim_local_nopad),
	_dim_local_last(other._dim_local_last)
{
    _data = other._data;
    other._data = nullptr;
}

template<typename Scalar>
SlicedVectorGpu<Scalar>&
SlicedVectorGpu<Scalar>::
operator=(const SlicedVectorGpu& other)
{
    // selbst
    if (this == &other) return *this;
    // fremd
    if(_data) delete[] _data;
	_my_comm = other._my_comm;
	_my_rank = other._my_rank;
	_num_nodes = other._num_nodes;
    _dim_global = other._dim_global;
    _dim_local = other._dim_local;
    _dim_local_nopad = other._dim_local_nopad;
	_dim_local_last = other._dim_local_last;
    _data = new Scalar[_dim_local];
    for (size_t i = 0; i < _dim_local; i++)
        _data[i] = other._data[i];
    return *this;
}

template<typename Scalar>
SlicedVectorGpu<Scalar>&
SlicedVectorGpu<Scalar>::
operator=(SlicedVectorGpu&& other)
{
    // selbst
    if (this == &other) return *this;
    // fremd
	_my_comm = other._my_comm;
	_my_rank = other._my_rank;
	_num_nodes = other._num_nodes;
    _dim_global = other._dim_global;
    _dim_local = other._dim_local;
    _dim_local_nopad = other._dim_local_nopad;
	_dim_local_last = other._dim_local_last;
    _data = other._data;
    other._data = nullptr;
    return *this;
}

template<typename Scalar>
void SlicedVectorGpu<Scalar>::
set_global(size_t pos, const Scalar& val)
{
    int affected_rank = pos / _dim_local_nopad;
    if (_my_rank == affected_rank)
        _data[pos - affected_rank*_dim_local_nopad] = val;
}

template<typename Scalar>
Scalar SlicedVectorGpu<Scalar>::
get_global(size_t pos) const
{
    int affected_rank = pos / _dim_local_nopad;
    Scalar val;
    if (_my_rank == affected_rank)
    {
        val = _data[pos - affected_rank *_dim_local_nopad];
    }
    MPI_SCALL(MPI_Bcast(&val,1,ScalarTraits<Scalar>::mpi_type,
                        affected_rank,_my_comm));
    return val;
}

template<typename Scalar>
void SlicedVectorGpu<Scalar>::print_local_data(std::ostream& out) const
{
	for (size_t i = 0; i < _dim_local; i++)
		out << i << ":\t" << _data[i] << std::endl;
}

template<typename Scalar>
typename SlicedVectorGpu<Scalar>::RealType
SlicedVectorGpu<Scalar>::
l2norm2_impl() const
{
    RealType res(0), res_global;
    for(size_t i=0; i<_dim_local; i++)
        res += ScalarTraits<Scalar>::abs2(_data[i]);
    MPI_SCALL(MPI_Allreduce(&res, &res_global, 1,
                            ScalarTraits<RealType>::mpi_type, MPI_SUM, _my_comm));
    return res_global;
}

template<typename Scalar>
typename SlicedVectorGpu<Scalar>::RealType
SlicedVectorGpu<Scalar>::
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

template<typename Scalar>
Scalar SlicedVectorGpu<Scalar>::
scal_prod_impl(const SlicedVectorGpu<Scalar>& other) const
{
    assert(_dim_global == other._dim_global);

    Scalar res(0), res_global;
    for(size_t i=0; i<_dim_local; i++)
        res += ScalarTraits<Scalar>::smult(_data[i], other._data[i]);
    MPI_SCALL(MPI_Allreduce(&res, &res_global, 1,
                            ScalarTraits<Scalar>::mpi_type, MPI_SUM, _my_comm));
    return res_global;
}

template<typename Scalar>
void SlicedVectorGpu<Scalar>::
axpy_impl(const Scalar& alpha, const SlicedVectorGpu<Scalar>& y)
{
    assert(_dim_global == y._dim_global);

    for(size_t i=0; i<_dim_local; i++)
        _data[i] = _data[i] + alpha*y._data[i];
}

template<typename Scalar>
void SlicedVectorGpu<Scalar>::
scal_impl(const Scalar& alpha)
{
   for(size_t i=0; i<_dim_local; i++)
        _data[i] *= alpha;
}


template<typename Scalar>
void SlicedVectorGpu<Scalar>::
swap_impl(SlicedVectorGpu &other)
{
    assert(_dim_global == other._dim_global);

    std::swap(_data,other._data);
}

template<typename Scalar>
void SlicedVectorGpu<Scalar>::
copy_impl(const SlicedVectorGpu &other)
{
    assert(_dim_global == other._dim_global);

    for(size_t i=0; i<_dim_local; i++) _data[i] = other._data[i];
}


/********** Spezialisierung der VectorTraits **********/


template<typename Scalar>
struct VectorTraits<SlicedVectorGpu<Scalar>>
{
    typedef typename ScalarTraits<Scalar>::RealType RealType;
    typedef Scalar ScalarType;
};

}

#endif // __SLICEDVECTORGPU_TPP_
