//#include <vector>
//#include "quadratur.hpp"
//#include "basis.hpp"

//#include "testfunctions.hpp"
#include "include/assemblefem.hpp"

namespace Icarus
{

double assembleFem::assemblyRHSLoad(std::vector<int>& e, std::vector<int>& A, mathfunction f)
{
    int n = e.size();
    double RHS(0.0);

    //TODO TOCHECK changed 02-24-16
    std::vector<double> X(27), Y(27), Z(27);
    std::vector<double> transx(27);
    std::vector<double> transy(27);
    std::vector<double> transz(27);
    transformation(_ax, h, transx);
    transformation(_ay, h, transy);
    transformation(_az, h, transz);
    //TODO TOCHECK changed 02-24-16

    for(int i = 0; i < n; i++)
    {
        //getQuadrature(e[i], "Name") = X, Y, Z, weigth;
        //TODO TOCHECK changed 02-24-16
        X = get_quadrature_xpoints(e[i], h, _ax, transx);
        Y = get_quadrature_xpoints(e[i], h, _ay, transy);
        Z = get_quadrature_xpoints(e[i], h, _az, transz);
        //TODO TOCHECK changed 02-24-16
        //get_quadrature_xpoints(e[i], X, h);
        int nqp = X.size();

        for(int q = 0; q<nqp; q++)
        {
            RHS += evaluate_Basis3d(e[i], A[i], X[q], Y[q], Z[q]) * f.eval(X[q], Y[q], Z[q]) * _weight[q];
        }
    }
    return RHS;
}

}//namespace Icarus
