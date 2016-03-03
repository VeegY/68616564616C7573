#ifndef _BASIS_HPP_
#define _BASIS_HPP_

#include <iostream>
#include <cmath>
#include <vector>

namespace Icarus
{
    double evaluate_Basis3d(size_t e, int A, double h, size_t Nx, size_t Ny, double X, double Y, double Z);

    std::vector<double> evaluate_gradient_Basis3d(size_t e, int A, double h, size_t Nx, size_t Ny, double X, double Y, double Z);

    double evaluate_Basis2d(size_t e, int A, double h, size_t Nx, size_t Ny, int type, double R1, double R2);
}

#endif
