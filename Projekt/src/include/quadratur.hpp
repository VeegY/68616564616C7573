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

#include <cmath>
#include <iostream>
#include <vector>
#include <math.h>


namespace Icarus
{
std::vector<double> get_weight(double c, double d);
void transformation(std::vector<double>& ai, double h, std::vector<double>& trans);
std::vector<double> get_quadrature_xpoints(int e, double h, std::vector<double>& ax, std::vector<double>& trans);
std::vector<double> get_quadrature_ypoints(int e, double h, std::vector<double>& ay, std::vector<double>& trans);
std::vector<double> get_quadrature_zpoints(int e, double h, std::vector<double>& az, std::vector<double>& trans);
}//namespace Icarus

#endif // __QUADRATUR_HPP_
