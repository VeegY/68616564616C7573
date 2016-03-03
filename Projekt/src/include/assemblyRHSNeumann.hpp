#ifndef __ASSEMBLYRHSNEUMANN_HPP_
#define __ASSEMBLYRHSNEUMANN_HPP_

#include <vector>
#include "basis.hpp"

#include "testfunctions.hpp"

namespace Icarus
{

double assemblyRHSNeumann(std::vector<size_t>& e, std::vector<int>& A, double h, size_t Nx, size_t Ny, int Ebene, math_function g=math_function(0))
{
    int n = e.size();
    double RHS(0.0);

    //TODO TOCHECK changed 02-24-16
    std::vector<double> X(27), Y(27), Z(27), weight(27);
    double aa(sqrt(0.6));
    std::vector<double> ax{-aa, -aa, -aa, -aa, -aa, -aa, -aa, -aa, -aa, 0, 0, 0, 0, 0, 0, 0, 0, 0, aa, aa, aa, aa, aa, aa, aa, aa, aa};  //x-Koordinaten der Gauss-Quadraturpunkte auf [-1,1]
    std::vector<double> ay{-aa, -aa, -aa, 0.5, 0, 0, aa, aa, aa, -aa, -aa, -aa, 0, 0, 0, aa, aa, aa, -aa, -aa, -aa, 0, 0, 0.25, aa, aa, aa}; //y-Koordinaten der Gauss-Quadraturpunkte auf [-1,1]
    std::vector<double> az{-aa, 0, aa, -aa, 0, aa, -aa, 0, aa, -aa, 0, aa, -aa, 0, aa, -aa, 0, aa, -aa, 0, aa, -aa, 0, aa, -aa, 0, aa};   //z-Koordinaten der Gauss-Quadraturpunkte auf [-1,1]
    std::vector<double> trans(27);
    //TODO TOCHECK changed 02-24-16

    for(int i = 0; i < n; i++)
    {
        //getQuadrature(e[i], "Name") = X, Y, Z, weigth;
        //TODO TOCHECK changed 02-24-16
        X = get_quadrature_xpoints(e[i], h, Nx, Ny, ax, trans);
        Y = get_quadrature_xpoints(e[i], h, Nx, Ny, ay, trans);
        Z = get_quadrature_xpoints(e[i], h, Nx, Ny, az, trans);
        //TODO TOCHECK changed 02-24-16
        //getQuadrature(e[i], "Name") = X, Y, Z, weight;

        int nqp = X.size();
            for(int q = 0; q<nqp; q++){
                //X-Y-Ebene
                if(Ebene == 1)
                {
                    RHS += evaluate_Basis2d(e[i], A[i], h, Nx, Ny Ebene, X[q], Y[q]) * g.eval(X[q], Y[q], Z[q]) * weight[q];
                }
                //X-Z-Ebene
                if(Ebene == 2)
                {
                    RHS += evaluate_Basis2d(e[i], A[i], h, Nx, Ny, Ebene, X[q], Z[q]) * g.eval(X[q], Y[q], Z[q]) * weight[q];
                }
                //Y-Z-Ebene
                if(Ebene == 3)
                {
                    RHS += evaluate_Basis2d(e[i], A[i], h, Nx, Ny, Ebene, Y[q], Z[q]) * g.eval(X[q], Y[q], Z[q]) * weight[q];
                }
            }
    }
    return RHS;
}

}//namespace Icarus

#endif
