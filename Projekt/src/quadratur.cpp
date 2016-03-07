#include "include/assemblefem.hpp"

namespace Icarus
{

std::vector<double> assembleFem::get_weight()  //Gibt die Gewichte der Quadratur als Vektor aus
{
    double c(5.0/9.0);
    double d(8.0/9.0);
    double e(c*c*c);
    double f((c*c)*d);
    double g(c*(d*d));
    return {e, f, e, f, g, f, e, f, e, f, g, f, g, d*d*d, g, f, g, f, e, f, e, f, g, f, e, f, e};
}

std::vector<double> assembleFem::get_weight_2d()  //Gibt die Gewichte der Quadratur als Vektor aus
{
    //TODO
    return {0.0};
}

std::vector<double> assembleFem::get_quadrature_xpoints(int e) //Berechnet die x-Koordinaten der Gauss-Quadraturpunkte für das Intervall für d Wuerfel mit Kantenlaenge hx*hy*hz
{
    //Berechnet die Translation und die Skalierung
    double e_x=getx(e);
    double a1=(1.0-sqrt(3.0/5.0))*h/2.0 + e_x;
    double a2=h/2.0 + e_x;
    double a3=(1.0+sqrt(3.0/5.0))*h/2.0 + e_x;

    //Bastelt den Vektor zusammen
    return std::vector<double>{a1,a1,a1,a1,a1,a1,a1,a1,a1, a2,a2,a2,a2,a2,a2,a2,a2,a2, a3,a3,a3,a3,a3,a3,a3,a3,a3};
}

std::vector<double> assembleFem::get_quadrature_ypoints(int e) //Berechnet die y-Koordinaten der Gauss-Quadraturpunkte fuer das Intervall für de Wuerfel mit Kantenlaenge hx*hy*h
{
    //Berechnet die Translation und die Skalierung
    double e_y=gety(e);
    double a1=(1.0-sqrt(3.0/5.0))*h/2.0 + e_y;
    double a2=h/2.0 + e_y;
    double a3=(1.0+sqrt(3.0/5.0))*h/2.0 + e_y;

    //Bastelt den Vektor zusammen
    return std::vector<double>{a1,a1,a1,a2,a2,a2,a3,a3,a3, a1,a1,a1,a2,a2,a2,a3,a3,a3, a1,a1,a1,a2,a2,a2,a3,a3,a3}; 
}


std::vector<double> assembleFem::get_quadrature_zpoints(int e) //Berechnet die z-Koordinaten der Gauss-Quadraturpunkte fuer das Intervall für den uerfel mit Kantenlaenge hx*hy*hz
{
    //Berechnet die Translation und die Skalierung
    double e_z=getz(e);
    double a1=(1.0-sqrt(3.0/5.0))*h/2.0 + e_z;
    double a2=h/2.0 + e_z;
    double a3=(1.0+sqrt(3.0/5.0))*h/2.0 + e_z;

    //Bastelt den Vektor zusammen
    return std::vector<double>{a1,a2,a3,a1,a2,a3,a1,a2,a3, a1,a2,a3,a1,a2,a3,a1,a2,a3, a1,a2,a3,a1,a2,a3,a1,a2,a3};
}

std::vector<double> assembleFem::get_quadrature_xpoints_2d(int e) //Berechnet die z-Koordinaten der Gauss-Quadraturpunkte fuer das Intervall für den uerfel mit Kantenlaenge hx*hy*hz
{
    //TODO
    return std::vector<double>{0.0};
}

std::vector<double> assembleFem::get_quadrature_ypoints_2d(int e) //Berechnet die z-Koordinaten der Gauss-Quadraturpunkte fuer das Intervall für den uerfel mit Kantenlaenge hx*hy*hz
{
    //TODO
    return std::vector<double>{0.0};
}

std::vector<double> assembleFem::get_quadrature_zpoints_2d(int e) //Berechnet die z-Koordinaten der Gauss-Quadraturpunkte fuer das Intervall für den uerfel mit Kantenlaenge hx*hy*hz
{
    //TODO
    return std::vector<double>{0.0};
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

