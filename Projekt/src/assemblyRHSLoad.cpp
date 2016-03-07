#include "include/assemblefem.hpp"

namespace Icarus
{

double assembleFem::assemblyRHSLoad(std::vector<int>& e, std::vector<int>& A, mathfunction f)
{
    int n = e.size();
    double RHS(0.0);

    std::vector<double> X(27), Y(27), Z(27);

    for(int i = 0; i < n; i++)
    {
        X = get_quadrature_xpoints(e[i]);
        Y = get_quadrature_ypoints(e[i]);
        Z = get_quadrature_zpoints(e[i]);
        int nqp = X.size();

        for(int q = 0; q<nqp; q++)
            RHS += evaluate_Basis3d(e[i], A[i], X[q], Y[q], Z[q]) * f.eval(X[q], Y[q], Z[q]) * _weight[q];
    }
    return RHS;
}

}//namespace Icarus
