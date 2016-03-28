#ifndef __DISTELLPACKMATRIXGPU_TPP_
#define __DISTELLPACKMATRIXGPU_TPP_

// nur für intellisense
#include "distellpackmatrixgpu.hpp"

template<typename mtype, typename vtype, typename rtype>
void gpu_ax_(mtype *data, const vtype* fvec, rtype* result, size_t *indices, size_t max_row_length, size_t dim_local);

template <typename Scalar>
void cleanupgpu(Scalar *data);

template<typename Scalar>
void alloc_unified(Scalar **fvec, size_t dim_fvec);

namespace Icarus
{

/**********  DistEllpackMatrixGpu  **********/


template <typename Scalar>
DistEllpackMatrixGpu<Scalar>::DistEllpackMatrixGpu(size_t dim_global, MPI_Comm my_comm) :
	_my_comm(my_comm),
    _dim_global(dim_global),
    _dim_local(0),
    _dim_local_nopad(0),
    _max_row_length(0),
    _indices(nullptr),
    _data(nullptr),
    _row_ptr(0),
    _col_ptr(0),
    _filled(false)
{
	// hole informationen über die mpi umgebung
	MPI_SCALL(MPI_Comm_size(_my_comm, &_num_nodes));
	MPI_SCALL(MPI_Comm_rank(_my_comm, &_my_rank));

    // geht die division genau auf
    if(_dim_global % _num_nodes == 0)
    {
        _dim_local = _dim_local_nopad = _dim_global/_num_nodes;
    }
    // funktioniert nur, wenn dim_global >> _num_nodes
    else
    {
        _dim_local = _dim_local_nopad = _dim_global/_num_nodes + 1;
        if(_my_rank == _num_nodes - 1)
            _dim_local = _dim_global - (_num_nodes - 1)*_dim_local_nopad;
        assert(_dim_local >= 0);
    }
}

template <typename Scalar>
DistEllpackMatrixGpu<Scalar>::
DistEllpackMatrixGpu(const DistEllpackMatrixGpu& other) :
_my_comm(other._my_comm),
_my_rank(other._my_rank),
_num_nodes(other._num_nodes),
_dim_global(other._dim_global),
_dim_local(other._dim_local),
_dim_local_nopad(other._dim_local_nopad),
_max_row_length(other._max_row_length),
_indices(nullptr),
_data(nullptr)
{
	try
	{
        alloc_unified(& _data, _max_row_length*_dim_local);
        alloc_unified(& _indices, _max_row_length*_dim_local);
	}
	catch (...)
	{
		LOG_ERROR("Memory allocation in DistEllpackMatrixGpu failed.");
	}
	for (size_t i = 0; i < _dim_local * _max_row_length; i++)
	{
		_indices[i] = other._indices[i];
		_data[i] = other._data[i];
	}
}

template <typename Scalar>
DistEllpackMatrixGpu<Scalar>::
DistEllpackMatrixGpu(DistEllpackMatrixGpu&& other) :
_my_comm(other._my_comm),
_my_rank(other._my_rank),
_num_nodes(other._num_nodes),
_dim_global(other._dim_global),
_dim_local(other._dim_local),
_dim_local_nopad(other._dim_local_nopad),
_max_row_length(other._max_row_length),
_indices(other._indices),
_data(other._data)
{
	other._indices = nullptr;
	other._data = nullptr;
}

template <typename Scalar>
DistEllpackMatrixGpu<Scalar>&
DistEllpackMatrixGpu<Scalar>::operator=(DistEllpackMatrixGpu&& other)
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
	_max_row_length = other._max_row_length;
	_indices = other._indices;
	_data = other._data;

	other._indices = nullptr;
	other._data = nullptr;

	return *this;
}

template <typename Scalar>
DistEllpackMatrixGpu<Scalar>&
DistEllpackMatrixGpu<Scalar>::operator=(const DistEllpackMatrixGpu& other)
{
	// selbst
	if (this == &other) return *this;

	// fremd
	if (_indices) cleanupgpu(_indices);
	if (_data) cleanupgpu(_data);

	_my_comm = other._my_comm;
	_my_rank = other._my_rank;
	_num_nodes = other._num_nodes;
	_dim_global = other._dim_global;
	_dim_local = other._dim_local;
	_dim_local_nopad = other._dim_local_nopad;
	_max_row_length = other._max_row_length;

	try
	{
	    alloc_unified(& _data, _max_row_length*_dim_local);
        alloc_unified(& _indices, _max_row_length*_dim_local);
	}
	catch (...)
	{
		LOG_ERROR("Memory allocation in DistEllpackMatrixGpu failed.");
	}
	for (size_t i = 0; i < _dim_local * _max_row_length; i++)
	{
		_indices[i] = other._indices[i];
		_data[i] = other._data[i];
	}

	return *this;
}


template <typename Scalar>
DistEllpackMatrixGpu<Scalar>::~DistEllpackMatrixGpu()
{
    if(_data) cleanupgpu(_data);
    if(_indices) cleanupgpu(_indices);
}

template <typename Scalar>
void DistEllpackMatrixGpu<Scalar>::prepare_sequential_fill(size_t max_row_length)
{
    // make sure the allocation has not yet taken place
    assert(!_data && !_indices);

    _max_row_length = max_row_length;
    try
    {
        alloc_unified(& _data, _max_row_length*_dim_local);
        alloc_unified(& _indices, _max_row_length*_dim_local);
    }
    catch(...)
    {
        LOG_ERROR("DistEllpackMatrixGpu: Memory allocation failed.");
    }
}

template <typename Scalar>
void DistEllpackMatrixGpu<Scalar>::sequential_fill(size_t colind, const Scalar& val)
{
    assert(!_filled);
    assert(_col_ptr < _max_row_length);
	assert(colind < _dim_global);

    _data[_col_ptr * _dim_local + _row_ptr] = val;
    _indices[_col_ptr * _dim_local + _row_ptr] = colind;
    _col_ptr++;
}

template <typename Scalar>
void DistEllpackMatrixGpu<Scalar>::end_of_row()
{
    assert(!_filled);

    // padding (maximally inefficient access pattern)
    for(; _col_ptr < _max_row_length; _col_ptr++)
    {
        _data[_col_ptr * _dim_local + _row_ptr] = PAD;
        _indices[_col_ptr * _dim_local + _row_ptr] = PAD;
    }

    // start new row
    _col_ptr = 0;
    if(++_row_ptr == _dim_local)
        _filled = true;
}



template <typename Scalar>
void DistEllpackMatrixGpu<Scalar>::mult_vec_impl(const VectorType& vec, VectorType& result) const
{
    assert(_dim_global == vec.get_dim_global());
    assert(_dim_global == result.get_dim_global());

    FullVectorGpu<Scalar> fvec(vec);

    gpu_ax_(_data, fvec.getDataPointer(), result.getDataPointer(), _indices, _max_row_length, _dim_local);
}


template <typename Scalar>
DistEllpackMatrixGpu<Scalar>
DistEllpackMatrixGpu<Scalar>::import_csr_file(const std::string &filename, MPI_Comm new_comm)
{
	// mpi umgebung des ziels
	int num_nodes, my_rank;
	MPI_SCALL(MPI_Comm_size(new_comm, &num_nodes));
	MPI_SCALL(MPI_Comm_rank(new_comm, &my_rank));

    std::ifstream file_colind, file_rowptr, file_vals;

    // drei teile der matrix öffnen
    file_colind.open(filename + "_colind.csr");
    file_rowptr.open(filename + "_rowptr.csr");
    file_vals.open(filename + "_vals.csr");
    if(!file_colind.is_open() || !file_rowptr.is_open() || !file_vals.is_open())
        LOG_ERROR("Failed to open ", filename, "_*.csr.");

    // bestimme globale dimension und #values
    size_t dim_global = get_num_lines(file_rowptr);
    size_t nnz = get_num_lines(file_vals);
    if(get_num_lines(file_colind) != nnz)
        LOG_ERROR("Inconsistent CSR data in ", filename, "_*.csr");

    DistEllpackMatrixGpu<Scalar> mat(dim_global, new_comm);

    // springe zur ersten zeile, die eingelesen wird
    go_to_line(file_rowptr,mat.first_row_on_node());

    const bool last_node = (my_rank == num_nodes - 1);

    // bestimme die maximale zeilenlänge auf der node
    size_t prev = 0, cur = 0, max_row_length = 0;
    file_rowptr >> prev;
    for(size_t i=0; i < mat.get_dim_local(); i++)
    {
        // global letzte zeile hat keinen nachfolger
        if(last_node && i == mat.get_dim_local()-1)
            cur = nnz;
        else
            file_rowptr >> cur;
        size_t row_len = cur-prev;
        if(row_len > max_row_length) max_row_length = row_len;
        prev = cur;
    }

    // bereite füllung auf mat-seite vor
    mat.prepare_sequential_fill(max_row_length);
    go_to_line(file_rowptr,mat.first_row_on_node());

    // vals und colind vorspulen
    file_rowptr >> prev;
    go_to_line(file_vals, prev);
    go_to_line(file_colind, prev);

    // zweiter durchgang: eigentliche füllung
    for(size_t i=0; i < mat.get_dim_local(); i++)
    {
        // global letzte zeile hat keinen nachfolger
        if(last_node && i == mat.get_dim_local()-1)
            cur = nnz;
        else
            file_rowptr >> cur;
        size_t row_len = cur-prev;

        for(size_t j=0; j < row_len; j++)
        {
            size_t colind;
            Scalar val;
            file_colind >> colind;
            file_vals >> val;
            mat.sequential_fill(colind,val);
        }
        mat.end_of_row();

        prev = cur;
    }

    if(!mat.is_filled())
        LOG_ERROR("csr_import failed unexpectedly.");

    LOG_INFO("Successfully loaded sparse matrix with dim=",dim_global,", nnz=",nnz,".");

    return mat;
}

template <typename Scalar>
DistEllpackMatrixGpu<Scalar>
DistEllpackMatrixGpu<Scalar>::precond_equi() const
{
    DistEllpackMatrixGpu<Scalar> Kinv(_dim_global);
    Scalar val;

    Kinv.prepare_sequential_fill(1);

    for(size_t i=0; i<_dim_local; i++)
    {
        val = 0;
        for(size_t j=0; j<_max_row_length; j++)
            val += ScalarTraits<Scalar>::abs(_data[j*_dim_local + i]);
        Kinv.sequential_fill(i, 1./val);
        Kinv.end_of_row();
    }
    return Kinv;
}

template <typename Scalar>
DistEllpackMatrixGpu<Scalar>
DistEllpackMatrixGpu<Scalar>::precond_jacobi() const
{
    DistEllpackMatrixGpu<Scalar> Kinv(_dim_global);
    size_t fron = first_row_on_node();
    Scalar val;

    Kinv.prepare_sequential_fill(1);

    for(size_t i=0; i<_dim_local; i++)
    {
          val = 1;
        for(size_t j=0; j<_max_row_length; j++)
        {
            if(_indices[j*_dim_local + i] == fron+i && _data[j*_dim_local + i] != 0)
            {
                val = 1./_data[j*_dim_local + i];
                break;
            }
        }
        Kinv.sequential_fill(i,val);
        Kinv.end_of_row();
    }

    return Kinv;
}

template <typename Scalar>
void DistEllpackMatrixGpu<Scalar>::print_local_data(std::ostream& os) const
{
    for(size_t i=0; i<_dim_local; i++)
    {
		os << i << ":\t";
        for(size_t j=0; j <_max_row_length; j++)
        {
            size_t pos = j*_dim_local + i;
            os << _data[pos] << "(" << _indices[pos] << ")\t";
        }
        os << std::endl;
    }
}

/**********  Passende MatrixTraits  **********/


template<typename Scalar>
struct MatrixTraits<DistEllpackMatrixGpu<Scalar>>
{
    typedef typename ScalarTraits<Scalar>::RealType RealType;
    typedef Scalar ScalarType;
    typedef SlicedVectorGpu<Scalar> VectorType;
};

}

#endif // __DISTELLPACKMATRIXGPU_TPP_
