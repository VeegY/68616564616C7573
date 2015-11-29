#ifndef __MATRIX_HPP_
#define __MATRIX_HPP_

#include <limits>
#include "slicedvector.hpp"
#include "fullvector.hpp"

namespace Icarus
{

template<typename T>
struct MatrixTraits;

// das interface für alle matrizen
template <class Child>
class Matrix
{
public:
    typedef typename MatrixTraits<Child>::RealType RealType;
    typedef typename MatrixTraits<Child>::ScalarType ScalarType;
    typedef typename MatrixTraits<Child>::VectorType VectorType;

protected:
    Matrix() { }

public:
    Child& leaf() { return static_cast<Child&>(*this); }
    const Child& leaf() const { return static_cast<const Child&>(*this); }

    void mult_vec(const Vector<VectorType>& x, Vector<VectorType>& res) const { leaf().mult_vec_impl(x.leaf(),res.leaf()); }

    size_t get_dim() const {return leaf().get_dim_impl();}
};

}


#endif // __MATRIX_HPP_
