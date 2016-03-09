#include "include/assemblefem.hpp"

namespace Icarus
{

double assembleFem::assemblyRHSLoad(mathfunction f)
{
    int n(_e.size());
    double RHS(0.0);

    std::vector<double> X(27), Y(27), Z(27);

    for (int i(0); i < n; i++)
    {
        X = get_quadrature_xpoints(_e[i]);
        Y = get_quadrature_ypoints(_e[i]);
        Z = get_quadrature_zpoints(_e[i]);

        for (int q(0); q < 27; q++)
            RHS += evaluate_Basis3d(_e[i], _A[i], X[q], Y[q], Z[q]) * f.eval(X[q], Y[q], Z[q]) * _weight[q];
    }
    return RHS;
}

}//namespace Icarus
