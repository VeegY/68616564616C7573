#ifndef __BICGSTABSOLVER_HPP_
#define __BICGSTABSOLVER_HPP_

#include <mpi.h>
#include "solver.hpp"
#include "matrix.hpp"
#include "vector.hpp"

namespace Icarus
{

/// \brief  Löst das LGS Ax=b mithilfe der Methode der bikonjugierten Gradienten
/// \tparam MatrixType  Typ der Matrix A und ggf. der Vorkonditionierer K1 und K2.
template <typename MatrixType>
class BiCgStabSolver : public Solver<BiCgStabSolver<MatrixType>>
{
    friend class Solver<BiCgStabSolver<MatrixType>>;

public:
    /// \brief  Mit der Matrix (bzgl. Speicherverteilung etc.) verträgliche Vektortyp.
    ///         Die rechte Seite b muss diesen Typ besitzen, die Lösung besitzt ebenfalls
    ///         diesen Typ.
    typedef typename MatrixType::VectorType VectorType;

    /// \brief  Typ der Einträge der Matrix A.
    typedef typename MatrixType::ScalarType ScalarType;

    /// \brief  Typ, den z.B. die Norm eines Vektors von Elementen vom Typ ScalarType hat.
    typedef typename MatrixType::RealType RealType;

private:
    RealType _tol;
    const MatrixType &_A;
    const VectorType &_b;
    const MatrixType *_K1inv, *_K2inv;

    MPI_Comm _comm;

public:
    /// \brief  Anzahl der Iterationen, nach der abgebrochen wird, wenn die Toleranz nicht
    ///         erreicht werden kann.
    static const long long MAX_ITER = 10000000000L;

    /// \brief  Toleranz, die ohne explizite Angabe angenommen wird
	static const RealType DEFAULT_TOL;


    /// \brief  Konstruktor
    //
    /// \param  A       Matrix A in Ax=b
    /// \param  b       Rechte Seite b in Ax=b
    /// \param  tol     Residuumsnorm, bei der abgebrochen werden soll
    /// \param  K1inv   Linksvorkonditionierer (nullptr für keinen)
    /// \param  K2inv   Rechtsvorkonditionierer (nullptr für keinen)
    BiCgStabSolver(
            const MatrixType& A,
            const VectorType& b,
            RealType tol = DEFAULT_TOL,
            const MatrixType* K1inv = nullptr,
            const MatrixType* K2inv = nullptr);

private:
    void solve_impl(VectorType& x0);

};

template <typename MatrixType>
const typename BiCgStabSolver<MatrixType>::RealType BiCgStabSolver<MatrixType>::DEFAULT_TOL(1e-9);

}

#include "bicgstabsolver.tpp"

#endif // __BICGSTABSOLVER_HPP
