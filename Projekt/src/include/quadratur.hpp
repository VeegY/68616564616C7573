/*
* Projekt:                  Studienprojekt TM 2015/16
* Dateiname:                quadratur.hpp
* Erstellt:                 28.1.16
* Autor / Ansprechpartner:  Justus
*
* Kurzbeschreibung:
* - Definiert das Interface Quadratur, das alle Operationen vorstellt,
*   die eine Quadratur implementieren muss.
*/

#ifndef __QUADRATUR_HPP_
#define __QUADRATUR_HPP_
#include <iostream>
#include <cmath>


namespace Icarus
{


class Quadratur
{
private:
    double c(5.0/9.0);

    double d(8.0/9.0);

    double a(sqrt(0.6));

public:

    double get_quadrature_xpoints(int e, double[] x_global, double h);

    double get_quadrature_ypoints(int e, double[] y_global, double h);

    double get_quadrature_zpoints(int e, double[] z_global, double h);

    double get_weight();
};

}

#endif // __QUADRATUR_HPP_
