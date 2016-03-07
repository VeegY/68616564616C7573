//#include <vector>
//#include "basis.hpp"

//#include "testfunctions.hpp"
#include "include/assemblefem.hpp"

namespace Icarus
{

double assembleFem::assemblyRHSNeumann(std::vector<int>& e, std::vector<int>& A, int Ebene, mathfunction g)
{
    int n = e.size();
    double RHS(0.0);

    std::vector<double> X(27), Y(27), Z(27);
    std::vector<double> trans(27);

    for(int i = 0; i < n; i++)
    {
        int nqp = X.size();
            for(int q = 0; q<nqp; q++){
                //X-Y-Ebene
                if(Ebene == 1)
                {
                    X = get_quadrature_xpoints_2d(e[i]);
                    Y = get_quadrature_ypoints_2d(e[i]);
                    RHS += evaluate_Basis2d(e[i], A[i], Ebene, X[q], Y[q]) * g.eval(X[q], Y[q], Z[q]) * _weight_2d[q];
                }
                //X-Z-Ebene
                if(Ebene == 2)
                {
                    X = get_quadrature_xpoints_2d(e[i]);
                    Z = get_quadrature_zpoints_2d(e[i]);
                    RHS += evaluate_Basis2d(e[i], A[i], Ebene, X[q], Z[q]) * g.eval(X[q], Y[q], Z[q]) * _weight_2d[q];
                }
                //Y-Z-Ebene
                if(Ebene == 3)
                {
                    Y = get_quadrature_ypoints_2d(e[i]);
                    Z = get_quadrature_zpoints_2d(e[i]);
                    RHS += evaluate_Basis2d(e[i], A[i], Ebene, Y[q], Z[q]) * g.eval(X[q], Y[q], Z[q]) * _weight_2d[q];
                }
            }
    }
    return RHS;
}

}//namespace Icarus
