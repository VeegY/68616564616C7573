#ifndef __SOLVER_HPP_
#define __SOLVER_HPP_

namespace Icarus
{

#include "vector.hpp"

template<typename T>
struct SolverTraits;

template <class Child>
class Solver
{
public:
    typedef typename SolverTraits<Child>::VectorType VectorType;

protected:
    Solver() { }

public:
    Child& leaf() { return static_cast<Child&>(*this); }
    const Child& leaf() const { return static_cast<const Child&>(*this); }

    void solve(Vector<VectorType>& dest) { leaf().solve_impl(dest.leaf()); }

};

}

#endif // __SOLVER_HPP_
