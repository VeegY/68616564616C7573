#ifndef __ASSEMBLYRHSLOAD_HPP_
#define __ASSEMBLYRHSLOAD_HPP_

#include <vector>
#include "quadratur.hpp"
#include "basis.hpp"

#include "testfunctions.hpp"

namespace Icarus
{

double assemblyRHSLoad(std::vector<int>& e, std::vector<int>& A, math_function f=math_function(0))
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
        X = get_quadrature_xpoints(e[i], h, ax, trans);
        Y = get_quadrature_xpoints(e[i], h, ay, trans);
        Z = get_quadrature_xpoints(e[i], h, az, trans);
        //TODO TOCHECK changed 02-24-16
        //get_quadrature_xpoints(e[i], X, h);
        int nqp = X.size();

        for(int q = 0; q<nqp; q++)
        {
            RHS += evaluate_Basis3d(e[i], A[i], X[q], Y[q], Z[q]) * f.eval(X[q], Y[q], Z[q]) * weight[q];
        }
    }
    return RHS;
}

}//namespace Icarus

#endif
