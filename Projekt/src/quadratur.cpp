#include "include/quadratur.hpp"
#include <cmath>
#include <iostream>
#include <vector>
#include <math.h>

namespace Icarus
{

std::vector<double> get_weight(double c, double d)
{
    return std::vector<double>{5.0};
}

void transformation(std::vector<double>& ai, double h, std::vector<double>& trans)
{
}

std::vector<double> get_quadrature_xpoints(int e, double h, std::vector<double>& ax, std::vector<double>& trans)
{
    return std::vector<double>{5.0};
}

std::vector<double> get_quadrature_ypoints(int e, double h, std::vector<double>& ay, std::vector<double>& trans)
{
    return std::vector<double>{5.0};
}

std::vector<double> get_quadrature_zpoints(int e, double h, std::vector<double>& az, std::vector<double>& trans)
{
    return std::vector<double>{5.0};
}

}//namespace Icarus
