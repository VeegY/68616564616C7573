#include "include/assemblefem.hpp"

namespace Icarus
{

//***** ***** ***** 3D ***** ***** *****//

std::vector<double> assembleFem::get_weight()  //Gibt die Gewichte der Quadratur als Vektor aus
{
    double c(5.0/9.0);
    double d(8.0/9.0);
    double e(c*c*c);
    double f((c*c)*d);
    double g(c*(d*d));
    return {e, f, e, f, g, f, e, f, e, f, g, f, g, d*d*d, g, f, g, f, e, f, e, f, g, f, e, f, e};
}

//Berechnet die x-Koordinaten der Gauss-Quadraturpunkte für das Referenzelement [0,h]^3
std::vector<double> assembleFem::get_quadrature_xpoints()
{
    //double a1=(1.0-sqrt(3.0/5.0))*_h/2.0;
    double a1(0.11270166537925831148*_h);
    double a2(_h/2.0);
    //double a3=(1.0+sqrt(3.0/5.0))*_h/2.0;
    double a3(0.88729833462074168852*_h);

    //Bastelt den Vektor zusammen
    return std::vector<double>{a1,a1,a1,a1,a1,a1,a1,a1,a1, a2,a2,a2,a2,a2,a2,a2,a2,a2, a3,a3,a3,a3,a3,a3,a3,a3,a3};
}

//Berechnet die y-Koordinaten der Gauss-Quadraturpunkte für das Referenzelement [0,h]^3
std::vector<double> assembleFem::get_quadrature_ypoints()
{
    //double a1=(1.0-sqrt(3.0/5.0))*_h/2.0;
    double a1(0.11270166537925831148*_h);
    double a2(_h/2.0);
    //double a3=(1.0+sqrt(3.0/5.0))*_h/2.0;
    double a3(0.88729833462074168852*_h);

    //Bastelt den Vektor zusammen
    return std::vector<double>{a1,a1,a1,a2,a2,a2,a3,a3,a3, a1,a1,a1,a2,a2,a2,a3,a3,a3, a1,a1,a1,a2,a2,a2,a3,a3,a3};
}

//Berechnet die z-Koordinaten der Gauss-Quadraturpunkte für das Referenzelement [0,h]^3
std::vector<double> assembleFem::get_quadrature_zpoints()
{
    //double a1=(1.0-sqrt(3.0/5.0))*_h/2.0;
    double a1(0.11270166537925831148*_h);
    double a2(_h/2.0);
    //double a3=(1.0+sqrt(3.0/5.0))*_h/2.0;
    double a3(0.88729833462074168852*_h);

    //Bastelt den Vektor zusammen
    return std::vector<double>{a1,a2,a3,a1,a2,a3,a1,a2,a3, a1,a2,a3,a1,a2,a3,a1,a2,a3, a1,a2,a3,a1,a2,a3,a1,a2,a3};
}

//***** ***** ***** 2D ***** ***** *****//

std::vector<double> assembleFem::get_weight_2d()  //Gibt die Gewichte der Quadratur als Vektor aus
{
    double c(5.0/9.0);
    double d(8.0/9.0);
    double e(c*c);
    double f(c*d);
    return {e,f,e, f,d*d,f, e,f,e};
}

std::vector<double> assembleFem::get_quadrature_xpoints_2d_1() //Berechnet die x-Koordinaten der Gauss-Quadraturpunkte für das Intervall für die Seitenflaeche (X-Y-Ebene)
{
    double a1((1.0-sqrt(3.0/5.0))*_h/2.0);
    double a2(_h/2.0);
    double a3((1.0+sqrt(3.0/5.0))*_h/2.0);

    //Bastelt den Vektor zusammen
    return std::vector<double>{a1,a1,a1, a2,a2,a2, a3,a3,a3};
}

std::vector<double> assembleFem::get_quadrature_ypoints_2d_1() //Berechnet die y-Koordinaten der Gauss-Quadraturpunkte für das Intervall für die Seitenflaeche (X-Y-Ebene)
{
    double a1((1.0-sqrt(3.0/5.0))*_h/2.0);
    double a2(_h/2.0);
    double a3((1.0+sqrt(3.0/5.0))*_h/2.0);

    //Bastelt den Vektor zusammen
    return std::vector<double>{a1,a2,a3, a1,a2,a3, a1,a2,a3};
}

std::vector<double> assembleFem::get_quadrature_xpoints_2d_2() //Berechnet die x-Koordinaten der Gauss-Quadraturpunkte für das Intervall für die Seitenflaeche (X-Z-Ebene)
{
    double a1((1.0-sqrt(3.0/5.0))*_h/2.0);
    double a2(_h/2.0);
    double a3((1.0+sqrt(3.0/5.0))*_h/2.0);

    //Bastelt den Vektor zusammen
    return std::vector<double>{a1,a1,a1, a2,a2,a2, a3,a3,a3};
}

std::vector<double> assembleFem::get_quadrature_zpoints_2d_2() //Berechnet die z-Koordinaten der Gauss-Quadraturpunkte für das Intervall für e Seitenflaeche (X-Z-Ebene)
{
    double a1((1.0-sqrt(3.0/5.0))*_h/2.0);
    double a2(_h/2.0);
    double a3((1.0+sqrt(3.0/5.0))*_h/2.0);

    //Bastelt den Vektor zusammen
    return std::vector<double>{a1,a2,a3, a1,a2,a3, a1,a2,a3};
}

std::vector<double> assembleFem::get_quadrature_ypoints_2d_3() //Berechnet die y-Koordinaten der Gauss-Quadraturpunkte für das Intervall für e Seitenflaeche (Y-Z-Ebene)
{
    double a1((1.0-sqrt(3.0/5.0))*_h/2.0);
    double a2(_h/2.0);
    double a3((1.0+sqrt(3.0/5.0))*_h/2.0);

    //Bastelt den Vektor zusammen
    return std::vector<double>{a1,a1,a1, a2,a2,a2, a3,a3,a3};
}

std::vector<double> assembleFem::get_quadrature_zpoints_2d_3() //Berechnet die z-Koordinaten der Gauss-Quadraturpunkte für das Intervall für d Seitenflaeche (Y-Z-Ebene)
{
    double a1=(1.0-sqrt(3.0/5.0))*_h/2.0;
    double a2=_h/2.0;
    double a3=(1.0+sqrt(3.0/5.0))*_h/2.0;

    //Bastelt den Vektor zusammen
    return std::vector<double>{a1,a2,a3, a1,a2,a3, a1,a2,a3};
}

}//namespace Icarus
