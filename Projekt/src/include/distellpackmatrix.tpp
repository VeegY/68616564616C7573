#ifndef __DISTMATRIXELLPACK_TPP_
#define __DISTMATRIXELLPACK_TPP_

// nur für intellisense
#include "distellpackmatrix.hpp"

namespace Icarus
{

/**********  DistEllpackMatrix  **********/


template <typename Scalar, int _num_nodes, int _first_node>
DistEllpackMatrix<Scalar, _num_nodes, _first_node>::DistEllpackMatrix(size_t dim_global) :
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
    // wenn weniger nodes vorhanden als angefordert, abbruch
    if(_last_node + 1 > MPI_HANDLER.get_n_procs())
        LOG_ERROR("DistEllpackMatrix is not compatible with node structure.");

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
}

template <typename Scalar, int _num_nodes, int _first_node>
DistEllpackMatrix<Scalar, _num_nodes, _first_node>::~DistEllpackMatrix()
{
    if(_data) delete[] _data;
    if(_indices) delete[] _indices;
}

template <typename Scalar, int _num_nodes, int _first_node>
void DistEllpackMatrix<Scalar, _num_nodes, _first_node>::prepare_sequential_fill(size_t max_row_length)
{
    // make sure the allocation has not yet taken place
    assert(!_data && !_indices);

    _max_row_length = max_row_length;
    try
    {
        _data = new Scalar[_dim_local * _max_row_length];
        _indices = new size_t[_dim_local * _max_row_length];
    }
    catch(...)
    {
        LOG_ERROR("DistEllpackMatrix: Memory allocation failed.");
    }
}

template <typename Scalar, int _num_nodes, int _first_node>
void DistEllpackMatrix<Scalar, _num_nodes, _first_node>::sequential_fill(size_t colind, const Scalar& val)
{
    assert(!_filled);
    assert(_col_ptr < _max_row_length);

    _data[_col_ptr * _dim_local + _row_ptr] = val;
    _indices[_col_ptr * _dim_local + _row_ptr] = colind;
    _col_ptr++;
}

template <typename Scalar, int _num_nodes, int _first_node>
void DistEllpackMatrix<Scalar, _num_nodes, _first_node>::end_of_row()
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


template <typename Scalar, int _num_nodes, int _first_node>
void DistEllpackMatrix<Scalar, _num_nodes, _first_node>::mult_vec_impl(const VectorType& vec, VectorType& result) const
{
    assert(_dim_global == vec.get_dim_global());
    assert(_dim_global == result.get_dim_global());

    // hole vec komplett in die node
    FullVector<Scalar> fvec(vec);

    // *************** BEGIN Durch CUDA-isierung ersetzen ************

    for(size_t row = 0; row < _dim_local; row++)
    {
        Scalar res(0);
        for(size_t col = 0; col < _max_row_length; col++)
        {
            size_t pos = col*_dim_local + row;
            if(_indices[pos] == PAD) continue;
            res += _data[pos] * fvec[_indices[pos]];
        }
        result.set_local(row, res);
    }

    // *************** END Durch CUDA-isierung ersetzen **************

}


/**********  Passende MatrixTraits  **********/


template<typename Scalar, int _num_nodes, int _first_node>
struct MatrixTraits<DistEllpackMatrix<Scalar, _num_nodes, _first_node>>
{
    typedef typename ScalarTraits<Scalar>::RealType RealType;
    typedef Scalar ScalarType;
    typedef SlicedVector<Scalar, _num_nodes, _first_node> VectorType;
};

}

#endif
