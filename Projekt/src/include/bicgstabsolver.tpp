#ifndef __BICGSTABSOLVER_TPP_
#define __BICGSTABSOLVER_TPP_

#include "bicgstabsolver.hpp"

namespace Icarus
{

template <typename MatrixType>
BiCgStabSolver<MatrixType>::BiCgStabSolver(
        const MatrixType& A,
        const VectorType& b,
        RealType tol,
        const MatrixType* K1inv,
        const MatrixType* K2inv) :
    _tol(tol),
    _A(A),
    _b(b),
    _K1inv(K1inv),
    _K2inv(K2inv),
    _comm(A.get_comm())
{
    assert(_A.get_dim() == _b.get_dim());
}

template <typename MatrixType>
void BiCgStabSolver<MatrixType>::solve_impl(VectorType& x0)
{
    assert(x0.get_dim() == _b.get_dim());
    const size_t dim = x0.get_dim();

    VectorType r_hat(dim, _comm), r(dim, _comm), nu(dim, _comm), s(dim, _comm), t(dim, _comm), p(dim, _comm);
    std::unique_ptr<VectorType> K1inv_t, K1inv_s, y, z;
    if(_K1inv)
    {
        K1inv_t.reset(new VectorType(dim, _comm));
        K1inv_s.reset(new VectorType(dim, _comm));
    }
    if(_K1inv || _K2inv)
    {
        y.reset(new VectorType(dim, _comm));
        z.reset(new VectorType(dim, _comm));
    }
    ScalarType rho, rho_, alpha, beta, omega;

    r.copy(_b);
    _A.mult_vec(x0,nu);
    r.axpy(-1,nu);

    r_hat.copy(r);

    rho_ = 1;
    alpha = 1;
    omega = 1;

    nu.clear();
    p.clear();

    for(long long i=0; i<MAX_ITER; i++)
    {
        rho = r_hat.scal_prod(r);
        beta = (rho/rho_)*(alpha/omega);

        p.scal(beta);
        p.axpy(1,r);
        p.axpy(-beta*omega,nu);

        if(_K1inv && _K2inv)
        {
            _K1inv->mult_vec(p,nu);
            _K2inv->mult_vec(nu,*y);
            _A.mult_vec(*y,nu);
        }
        else if(_K2inv)
        {
            _K2inv->mult_vec(p,*y);
            _A.mult_vec(*y,nu);
        }
        else if(_K1inv)
        {
            _K1inv->mult_vec(p,*y);
            _A.mult_vec(*y,nu);
        }
        else
            _A.mult_vec(p, nu);

        alpha = rho/r_hat.scal_prod(nu);

        s.copy(r);
        s.axpy(-alpha, nu);
        LOG_DEBUG("BiCgStab: After ",i+1," its, sq_tol = ", s.l2norm2());
        if(s.l2norm2() < _tol)
        {
            x0.axpy(alpha, p);
            return;
        }

        if(_K1inv && _K2inv)
        {
            _K1inv->mult_vec(s,t);
            _K2inv->mult_vec(t,*z);
            _A.mult_vec(*z,t);
        }
        else if(_K1inv)
        {
            _K1inv->mult_vec(s,*z);
            _A.mult_vec(*z,t);
        }
        else if(_K2inv)
        {
            _K2inv->mult_vec(s,*z);
            _A.mult_vec(*z,t);
        }
        else
            _A.mult_vec(s,t);

        if(_K1inv)
        {
            _K1inv->mult_vec(t,*K1inv_t);
            _K1inv->mult_vec(s,*K1inv_s);
            omega = K1inv_t->scal_prod(*K1inv_s) / K1inv_t->l2norm2();
        }
        else
            omega = t.scal_prod(s) / t.l2norm2();

        if(y) x0.axpy(alpha,*y);
        else x0.axpy(alpha,p);
        if(z) x0.axpy(omega,*z);
        else x0.axpy(omega,s);

        r.copy(s);
        r.axpy(-omega, t);

        // shift index
        rho_ = rho;
    }

    LOG_ERROR("BiCgStab: Not converged after ", MAX_ITER, " iterations.");
}


template <typename MatrixT>
struct SolverTraits<BiCgStabSolver<MatrixT>>
{
    typedef typename MatrixT::VectorType VectorType;
};

}//namespace Icarus

#endif // __BICGSTABSOLVER_TPP_
