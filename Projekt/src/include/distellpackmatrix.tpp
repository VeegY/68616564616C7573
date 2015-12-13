#ifndef __DISTELLPACKMATRIX_TPP_
#define __DISTELLPACKMATRIX_TPP_

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
            if(_data[pos] == PAD) continue;
            res += _data[pos] * fvec[_indices[pos]];
        }
        result.set_local(row, res);
    }

    // *************** END Durch CUDA-isierung ersetzen **************

}


template <typename Scalar, int _num_nodes, int _first_node>
DistEllpackMatrix<Scalar, _num_nodes, _first_node>
DistEllpackMatrix<Scalar, _num_nodes, _first_node>::import_csr_file(const std::string &filename)
{
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

    DistEllpackMatrix<Scalar, _num_nodes, _first_node> mat(dim_global);

    // wenn wir nicht am füllen beteiligt sind, fertig
    const int this_node = MPI_HANDLER.get_my_rank();
    if(this_node >= _first_node + _num_nodes)
        return mat;

    // springe zur ersten zeile, die eingelesen wird
    go_to_line(file_rowptr,mat.first_row_on_node());

    const bool last_node = (this_node == _first_node + _num_nodes - 1);

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


template <typename Scalar, int _num_nodes, int _first_node>
void DistEllpackMatrix<Scalar, _num_nodes, _first_node>::print_local_data(std::ostream& os) const
{
    for(size_t i=0; i<_dim_local; i++)
    {
        for(size_t j=0; j <_max_row_length; j++)
        {
            size_t pos = j*_dim_local + i;
            os << _data[pos] << "(" << _indices[pos] << ")\t";
        }
        os << std::endl;
    }
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

#endif // __DISTELLPACKMATRIX_TPP_
