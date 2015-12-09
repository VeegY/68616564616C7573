#ifndef __BICGSTABSOLVER_HPP_
#define __BICGSTABSOLVER_HPP_

#include "solver.hpp"
#include "matrix.hpp"
#include "vector.hpp"

namespace Icarus
{

template <typename MatrixType>
class BiCgStabSolver : public Solver<BiCgStabSolver<MatrixType>>
{
    friend class Solver<BiCgStabSolver<MatrixType>>;

public:
    typedef typename MatrixType::VectorType VectorType;
    typedef typename MatrixType::ScalarType ScalarType;
    typedef typename MatrixType::RealType RealType;

private:
    RealType _tol;
    const MatrixType &_A, *_K1inv, *_K2inv;
    const VectorType &_b;

public:
    static const int MAX_ITER = 1000;
    static constexpr RealType DEFAULT_TOL = RealType(1e-9);

    BiCgStabSolver(
            const MatrixType& A,
            const VectorType& b,
            RealType tol = DEFAULT_TOL,
            const MatrixType* K1inv = nullptr,
            const MatrixType* K2inv = nullptr);

private:
    void solve_impl(VectorType& x0);

};

}

#include "bicgstabsolver.tpp"

#endif // __BICGSTABSOLVER_HPP
