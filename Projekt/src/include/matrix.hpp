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
    Child& _leaf;
public:
    typedef typename MatrixTraits<Child>::RealType RealType;
    typedef typename MatrixTraits<Child>::ScalarType ScalarType;
    typedef typename MatrixTraits<Child>::VectorType VectorType;

protected:
    Matrix() : _leaf(*static_cast<Child*>(this)) { }

public:
    void mult_vec(const VectorType& x, VectorType& res) const { _leaf.mult_vec_impl(x,res); }

};

}


#endif // __MATRIX_HPP_
