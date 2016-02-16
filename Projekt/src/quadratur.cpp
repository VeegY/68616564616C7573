#include "include/quadratur.hpp"
#include <cmath>
#include <iostream>
#include <vector>
#include <math.h>

namespace Icarus
{

std::vector<double> get_weight(double c, double d)  //Gibt die Gewichte der Quadratur aus
{
    double e(c*c*c);
    double f((c*c)*d);
    double g(c*(d*d));
    return {e, f, e, f, g, f, e, f, e, f, g, f, g, d*d*d, g, f, g, f, e, f, e, f, g, f, e, f, e};
}

void transformation(std::vector<double>& ai, double h, std::vector<double>& trans) //Berechnet 1/2*h*ax[i]+1/2*h zur Transformation auf das Intervall der Länge h
{
    for(int i=0; i<27; ++i)
    {
        trans[i]=(ai[i]+1)*0.5*h;
    }
}

std::vector<double> get_quadrature_xpoints(int e, double h, std::vector<double>& ax, std::vector<double>& trans) //Berechnet die x-Koordinaten der Gauss-Quadraturpunkte für das Intervall für den Würfel mit Kantenlänge hx*hy*hz
{
    std::vector<double> x_global(27);
    //e_x=getXcoordinate(e);
    std::vector<double> x_local{0, 0, h, h, 0, 0, h, h, 0, 0.5*h, h, 0.5*h, 0, 0,
        h, h, 0, 0.5*h, h, 0.5*h, 0.5*h, 0, 0.5*h, h, 0.5*h, 0.5*h, 0.5*h};
    transformation(ax, h, trans);
    for(int j=0; j<27; ++j)
    {
        x_global[j]=trans[j]+1.0+x_local[j];
    }
    return x_global;
}

std::vector<double> get_quadrature_ypoints(int e, double h, std::vector<double>& ay, std::vector<double>& trans) //Berechnet die y-Koordinaten der Gauss-Quadraturpunkte für das Intervall für den Würfel mit Kantenlänge hx*hy*hz
{
    std::vector<double> y_global(27);
    //e_y=getYcoordinate(e);
    std::vector<double> y_local{h, 0, 0, h, h, 0, 0, h, 0.5*h, h, 0.5*h, h, h,
        0, 0, h, 0.5*h, 0, 0.5*h, h, 0.5*h, 0.5*h, 0, 0.5*h, h, 0.5*h, 0.5*h};
    transformation(ay, h, trans);
    for(int l=0; l<27; ++l)
    {
        y_global[l]=trans[l]+1.0+y_local[l];
    }
    return y_global;
}



std::vector<double> get_quadrature_zpoints(int e, double h, std::vector<double>& az, std::vector<double>& trans) //Berechnet die z-Koordinaten der Gauss-Quadraturpunkte für das Intervall für den Würfel mit Kantenlänge hx*hy*hz
{
    std::vector<double> z_global(27);
    //e_z=getZcoordinate(e);
    std::vector<double> z_local{0, 0, 0, 0, h, h, h, h, 0, 0, 0, 0,
        0.5*h, 0.5*h, 0.5*h, 0.5*h, h, h, h, h, 0,  0.5*h, 0.5*h, 0.5*h, 0.5*h, h,  0.5*h};
    transformation(az, h, trans);
    for(int r=0; r<27; ++r)
    {
        z_global[r]=trans[r]+1.0+z_local[r];
    }
    return z_global;
}

/*
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
*/

}//namespace Icarus
