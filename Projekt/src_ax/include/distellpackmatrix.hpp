#ifndef __DISTELLPACKMATRIX_HPP_
#define __DISTELLPACKMATRIX_HPP_

#include "matrix.hpp"
#include "utility.hpp"
#include <fstream>
#include <iterator>
#include <algorithm>

namespace Icarus
{

// quadratische matrix, global zeilenweise verteilt, lokal
// in ellpack gespeichert
template <typename Scalar, int _num_nodes, int _first_node = 0>
class DistEllpackMatrix : public Matrix<DistEllpackMatrix<Scalar, _num_nodes, _first_node>>
{
    friend class Matrix<DistEllpackMatrix<Scalar, _num_nodes, _first_node>>;

    // Mit PAD wird das padding durchgeführt
    static constexpr int PAD = 0;

    static constexpr int _last_node = _first_node + _num_nodes - 1;

    size_t _dim_global, _dim_local, _dim_local_nopad, _max_row_length;

    size_t * _indices;

    Scalar* _data;

    // hilfsvariablen zum sequentiellen füllen der matrix
    size_t _row_ptr, _col_ptr;
    bool _filled;

public:

    typedef typename MatrixTraits<DistEllpackMatrix<Scalar, _num_nodes, _first_node>>::VectorType VectorType;

    DistEllpackMatrix(size_t dim_global);

    ~DistEllpackMatrix();

    //TODO: copy, move, assignment, move assignment

    size_t get_dim_local() const { return _dim_local; }

    size_t get_dim_local_nopad() const { return _dim_local_nopad; }

    size_t get_dim_global() const { return _dim_global; }

    //TODO: get, set local/global

    void prepare_sequential_fill(size_t max_row_length);

    void sequential_fill(size_t colind, const Scalar& val);

    void end_of_row();

    bool is_filled() const { return _filled; }

    size_t first_row_on_node() const { return (MPI_HANDLER.get_my_rank() - _first_node) * _dim_local_nopad; }

    static DistEllpackMatrix import_csr_file(const std::string& filename);

    void print_local_data(std::ostream &os) const;

private:

    size_t get_dim_impl() const {return _dim_global;}

    void mult_vec_impl(const VectorType& vec, VectorType& result) const;
};

}

#include "distellpackmatrix.tpp"




#endif // __DISTELLPACKMATRIX_HPP_
