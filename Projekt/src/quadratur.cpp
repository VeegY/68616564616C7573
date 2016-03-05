//#include "include/quadratur.hpp"
//#include "include/getxyz.hpp"
#include <cmath>
#include <iostream>
#include <vector>
#include <math.h>

#include "include/assemblefem.hpp"

namespace Icarus
{

std::vector<double> assembleFem::get_weight(double c, double d)  //Gibt die Gewichte der Quadratur als Vektor aus
{
    double e(c*c*c);
    double f((c*c)*d);
    double g(c*(d*d));
    return {e, f, e, f, g, f, e, f, e, f, g, f, g, d*d*d, g, f, g, f, e, f, e, f, g, f, e, f, e};
}

void assembleFem::transformation(std::vector<double>& ai, std::vector<double>& trans) //Berechnet Transformation auf das Intervall der Länge h
{
    for(int i=0; i<27; ++i)
    {
        trans[i]=(ai[i]+1)*0.5*h;
    }
}

std::vector<double> assembleFem::get_quadrature_xpoints(int e, std::vector<double>& trans) //Berechnet die x-Koordinaten der Gauss-Quadraturpunkte für das Intervall für den Würfel mit Kantenlänge hx*hy*hz
{

    std::vector<double> x_global(27);
    int e_x=getx(e);
    std::vector<double> x_local{0, 0, h, h, 0, 0, h, h, 0, 0.5*h, h, 0.5*h, 0, 0,
        h, h, 0, 0.5*h, h, 0.5*h, 0.5*h, 0, 0.5*h, h, 0.5*h, 0.5*h, 0.5*h};
    transformation(_ax, trans); //Transformation auf das Intervall mit der Länge h
    for(int j=0; j<27; ++j)
    {
        x_global[j]=trans[j]+e_x+x_local[j]; 
    }
    return x_global;
}

std::vector<double> assembleFem::get_quadrature_ypoints(int e, std::vector<double>& trans) //Berechnet die y-Koordinaten der Gauss-Quadraturpunkte für das Intervall für den Würfel mit Kantenlänge hx*hy*hz
{
    std::vector<double> y_global(27);
    int e_y=gety(e);
    std::vector<double> y_local{h, 0, 0, h, h, 0, 0, h, 0.5*h, h, 0.5*h, h, h,
        0, 0, h, 0.5*h, 0, 0.5*h, h, 0.5*h, 0.5*h, 0, 0.5*h, h, 0.5*h, 0.5*h};
    transformation(_ay, trans); //Transformation auf das Intervall mit der Länge h
    for(int l=0; l<27; ++l)
    {
        y_global[l]=trans[l]+e_y+y_local[l];
    }
    return y_global;
}


std::vector<double> assembleFem::get_quadrature_zpoints(int e, std::vector<double>& trans) //Berechnet die z-Koordinaten der Gauss-Quadraturpunkte für das Intervall für den Würfel mit Kantenlänge hx*hy*hz
{
    std::vector<double> z_global(27);
    int e_z=getz(e);
    std::vector<double> z_local{0, 0, 0, 0, h, h, h, h, 0, 0, 0, 0,
        0.5*h, 0.5*h, 0.5*h, 0.5*h, h, h, h, h, 0,  0.5*h, 0.5*h, 0.5*h, 0.5*h, h,  0.5*h};
    transformation(_az, trans); //Transformation auf das Intervall mit der Länge h
    for(int r=0; r<27; ++r)
    {
        z_global[r]=trans[r]+e_z+z_local[r]; 
    }
    return z_global;
}



}//namespace Icarus

/*Beispielhafte main-Funktion zum Verständnis, die einzelnen Funktionen benutzt/initialisiert werden müssen:

int main()
{
    double c(5.0/9.0);
    double d(8.0/9.0);
    std::vector<double> weight{(Icarus::get_weight(c, d))}; //Vektor der Gewichte kann jetzt mit "weight" benutzt werden


    double h(0.01); // Beispielhafte Schrittweite

    double a(sqrt(0.6));
    std::vector<double> ax{-a, -a, -a, -a, -a, -a, -a, -a, -a, 0, 0, 0, 0, 0, 0, 0, 0, 0, a, a, a, a, a, a, a, a, a};  //x-Koordinaten der Gauss-Quadraturpunkte auf [-1,1]
    std::vector<double> ay{-a, -a, -a, 0.5, 0, 0, a, a, a, -a, -a, -a, 0, 0, 0, a, a, a, -a, -a, -a, 0, 0, 0.25, a, a, a}; //y-Koordinaten der Gauss-Quadraturpunkte auf [-1,1]
    std::vector<double> az{-a, 0, a, -a, 0, a, -a, 0, a, -a, 0, a, -a, 0, a, -a, 0, a, -a, 0, a, -a, 0, a, -a, 0, a};   //z-Koordinaten der Gauss-Quadraturpunkte auf [-1,1]


    int e(1); //Beispielhaftes e


    std::vector<double> ai(27);//Wird benutzt in der Funktion transformation
    std::vector<double> trans(27);//Wird bei der Funktion transformation als Vektor ausgegeben

    std::vector<double> x_global(27);//X-Koordinaten der entgültigen Quadraturpunkte
    Icarus::get_quadrature_xpoints(e, h, ax, trans);
    std::vector<double> y_global(27);//Y-Koordinaten der entgültigen Quadraturpunkte
    Icarus::get_quadrature_ypoints(e, h, ay, trans);
    std::vector<double> z_global(27); //Z-Koordinaten der entgültigen Quadraturpunkte
    Icarus::get_quadrature_zpoints(e, h, az, trans);


    return 0;
}




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

