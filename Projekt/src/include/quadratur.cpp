#ifndef __QUADRATUR_CPP_
#define __QUADRATUR_CPP_

#include "quadratur.hpp"
#include <cmath>
#include <iostream>
#include <vector>
#include <math.h>

namespace Icarus
{



std::vector<double> get_weight(double c, double d)
{
    return {c*c*c, (c*c)*d, c*c*c, (c*c)*d, c*(d*d), (c*c)*d, c*c*c,
    (c*c)*d, c*c*c, (c*c)*d, c*(d*d), (c*c)*d, c*(d*d), d*d*d,
    c*(d*d),(c*c)*d, c*(d*d), (c*c)*d, c*c*c,
    (c*c)*d, c*c*c, (c*c)*d, c*(d*d), (c*c)*d, c*c*c, (c*c)*d, c*c*c};
}


std::vector<double> ax={-a, -a, -a, -a, -a, -a, -a, -a, -a, 0, 0, 0, 0, 0, 0, 0, 0, 0, a, a, a, a, a, a, a, a, a};  //x-Koordinaten der Gauss-Quadraturpunkte auf [-1,1]

std::vector<double> ay={-a, -a, -a, 0.5, 0, 0, a, a, a, -a, -a, -a, 0, 0, 0, a, a, a, -a, -a, -a, 0, 0, 0.25, a, a, a}; //y-Koordinaten der Gauss-Quadraturpunkte auf [-1,1]

std::vector<double> az={-a, 0, a, -a, 0, a, -a, 0, a, -a, 0, a, -a, 0, a, -a, 0, a, -a, 0, a, -a, 0, a, -a, 0, a};   //z-Koordinaten der Gauss-Quadraturpunkte auf [-1,1]

std::vector<double> ai(27);
std::vector<double> trans(27);

void transformation(std::vector<double>& ai, double h, std::vector<double>& trans) //Berechnet 1/2*h*ax[i]+1/2*h zur Transformation auf das Intervall der Länge h
{
    for(int i=0; i<27; ++i)
    {
        trans[i]=ai[i]*0.5*h+0.5*h;
    }
}



std::vector<double> get_quadrature_xpoints(int e, double[] x_global, double h) //Berechnet die x-Koordinaten der Gauss-Quadraturpunkte für das Intervall für den Würfel mit Kantenlänge hx*hy*hz
{

    e_x=getXcoordinate(e);
    double x_local[27]={0, 0, h, h, 0, 0, h, h, 0, 0.5*h, h, 0.5*h, 0, 0, h, h, 0, 0.5*h, h, 0.5*h, 0.5*h, 0, 0.5*h, h, 0.5*h, 0.5*h, 0.5*h};

    for(int j=0; j<27; ++j)
    {
        x_global[j]=transformation(ax[j])+e_x+x_local[j];
    }
    return x_global;
}

double get_quadrature_ypoints(int e, double[] y_global, double h) //Berechnet die y-Koordinaten der Gauss-Quadraturpunkte für das Intervall für den Würfel mit Kantenlänge hx*hy*hz
{
    e_y=getYcoordinate(e);
    double y_local[27]={h, 0, 0, h, h, 0, 0, h, 0.5*h, h, 0.5*h, h, h, 0, 0, h, 0.5*h, 0, 0.5*h, h, 0.5*h, 0.5*h, 0, 0.5*h, h, 0.5*h, 0.5*h};
    for(int l=0; l<27; ++l)
    {
        y_global[l]=transformation(ay[l])+e_y+y_local[l];
    }
    return y_global;
}

double get_quadrature_zpoints(int e, double[] z_global, double h) //Berechnet die z-Koordinaten der Gauss-Quadraturpunkte für das Intervall für den Würfel mit Kantenlänge hx*hy*hz
{
    e_z=getZcoordinate(e);
    double z_local[27]={0, 0, 0, 0, h, h, h, h, 0, 0, 0, 0, 0.5*h, 0.5*h, 0.5*h, 0.5*h, h, h, h, h, 0,  0.5*h, 0.5*h, 0.5*h, 0.5*h, h,  0.5*h};
    for(int m=0; 0<27; ++m)
    {
        z_global[m]=transformation(az[m])+e_z+z_local[m];
    }
    return z_global;
}

}

#endif//__QUADRATUR_CPP_
