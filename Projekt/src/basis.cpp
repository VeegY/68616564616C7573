//TODO TODISCUSS:
// immer direkt return ...; ? also zwsp und x0,y0,z0 umgehen

#include "include/assemblefem.hpp"

namespace Icarus
{

double assembleFem::evaluate_Basis3d(int e, int A, double X, double Y, double Z)
{
    double x0(getx(e));
    double y0(gety(e));
    double z0(getz(e));
    double zwsp(0.0);

    switch(A)
    {
        case 0: zwsp = ((x0 + h - X) / h) * ((y0 + h - Y) / h) * ((z0 + h - Z) / h);
                break;
        case 1: zwsp = ((x0 - X) / -h) * ((y0 + h - Y) / h) * ((z0 + h - Z) / h);
                break;
        case 2: zwsp = ((x0 - X) / -h) * ((y0 - Y) / -h) * ((z0 + h - Z) / h);
                break;
        case 3: zwsp = ((x0 + h - X) / h) * ((y0 - Y) / -h) * ((z0 + h - Z) / h);
                break;
        case 4: zwsp = ((x0 + h - X) / h) * ((y0 + h - Y) / h) * ((z0 - Z) / -h);
                break;
        case 5: zwsp = ((x0 - X) / -h) * ((y0 + h - Y) / h) * ((z0 - Z) / -h);
                break;
        case 6: zwsp = ((x0 - X) / -h) * ((y0 - Y) / -h) * ((z0 - Z) / -h);
                break;
        case 7: zwsp = ((x0 + h - X) / h) * ((y0 - Y) / -h) * ((z0 - Z) / -h);
                break;
        default: std::cout << "Fehler: kein lokaler Knoten" << std::endl;
    }
    return zwsp;
}


std::vector<double> assembleFem::evaluate_gradient_Basis3d(int e, int A, double X, double Y, double Z)
{
    double x0(getx(e));
    double y0(gety(e));
    double z0(getz(e));
    std::vector<double> zwsp{0.0, 0.0, 0.0};

    switch(A)
    {
        case 0: zwsp[0] = (-1/h) * ((y0 + h - Y) / h) * ((z0 + h - Z) / h);
                zwsp[1] = (-1/h) * ((x0 + h - X) / h) * ((z0 + h - Z) / h);
                zwsp[2] = (-1/h) * ((x0 + h - X) / h) * ((y0 + h - Y) / h);
                break;
        case 1: zwsp[0] = (1/h) * ((y0 + h - Y) / h) * ((z0 + h - Z) / h);
                zwsp[1] = (-1/h) * ((x0 - X) / -h) * ((z0 + h - Z) / h);
                zwsp[2] = (-1/h) * ((x0 - X) / -h) * ((y0 + h - Y) / h);
                break;
        case 2: zwsp[0] = (1/h) * ((y0 - Y) / -h) * ((z0 + h - Z) / h);
                zwsp[1] = (1/h) * ((x0 - X) / -h) * ((z0 + h - Z) / h);
                zwsp[2] = (-1/h) * ((x0 - X) / -h) * ((y0 - Y) / -h);
                break;
        case 3: zwsp[0] = (-1/h) * ((y0 - Y) / -h) * ((z0 + h - Z) / h);
                zwsp[1] = (1/h) * ((x0 + h - X) / h) * ((z0 + h - Z) / h);
                zwsp[2] = (-1/h) * ((x0 + h - X) / h) * ((y0 - Y) / -h);
                break;
        case 4: zwsp[0] = (-1/h) * ((y0 + h - Y) / h) * ((z0 - Z) / -h);
                zwsp[1] = (-1/h) * ((x0 + h - X) / h) * ((z0 - Z) / -h);
                zwsp[2] = (1/h) * ((x0 + h - X) / h) * ((y0 + h - Y) / h);
                break;
        case 5: zwsp[0] = (1/h) * ((y0 + h - Y) / h) * ((z0 - Z) / -h);
                zwsp[1] = (-1/h) * ((x0 - X) / -h) * ((z0 - Z) / -h);
                zwsp[2] = (1/h) * ((x0 - X) / -h) * ((y0 + h - Y) / h);
                break;
        case 6: zwsp[0] = (1/h) * ((y0 - Y) / -h) * ((z0 - Z) / -h);
                zwsp[1] = (1/h) * ((x0 - X) / -h) * ((z0 - Z) / -h);
                zwsp[2] = (1/h) * ((x0 - X) / -h) * ((y0 - Y) / -h);
                break;
        case 7: zwsp[0] = (-1/h) * ((y0 - Y) / -h) * ((z0 - Z) / -h);
                zwsp[1] = (1/h) * ((x0 + h - X) / h) * ((z0 - Z) / -h);
                zwsp[2] = (1/h) * ((x0 + h - X) / h) * ((y0 - Y) / -h);
                break;
        default: std::cout << "Fehler: kein lokaler Knoten" << std::endl;
    }
    return zwsp;
}

double assembleFem::evaluate_Basis2d_1(int e, int A, double R1, double R2)
{
    double x0(getx(e));
    double y0(gety(e));
    double zwsp(0.0);

    switch (A)
    {
        case 0: zwsp = ((x0 + h - R1) / h) * ((y0 + h - R2) / h);
                break;
        case 1: zwsp = ((x0 - R1) / -h) * ((y0 + h - R2) / h);
                break;
        case 2: zwsp = ((x0 - R1) / -h) * ((y0 - R2) / -h);
                break;
        case 3: zwsp = ((x0 + h - R1) / h) * ((y0 - R2) / -h);
                break;
        default: std::cout << "Fehler: kein lokaler Knoten" << std::endl;
    }
    return zwsp;
}

double assembleFem::evaluate_Basis2d_2(int e, int A, double R1, double R2)
{
    double x0(getx(e));
    double z0(getz(e));
    double zwsp(0.0);

    switch (A)
    {
        case 0: zwsp = ((x0 + h - R1) / h) * ((z0 + h - R2) / h);
                break;
        case 1: zwsp = ((x0 - R1) / -h) * ((z0 + h - R2) / h);
                break;
        case 2: zwsp = ((x0 - R1) / -h) * ((z0 - R2) / -h);
                break;
        case 3: zwsp = ((x0 + h - R1) / h) * ((z0 - R2) / -h);
                break;
        default: std::cout << "Fehler: kein lokaler Knoten" << std::endl;
    }
    return zwsp;
}

double assembleFem::evaluate_Basis2d_3(int e, int A, double R1, double R2)
{
    double y0(gety(e));
    double z0(getz(e));
    double zwsp(0.0);

    switch (A)
    {
        case 0: zwsp = ((y0 + h - R1) / h) * ((z0 + h - R2) / h);
                break;
        case 1: zwsp = ((y0 - R1) / -h) * ((z0 + h - R2) / h);
                break;
        case 2: zwsp = ((y0 - R1) / -h) * ((z0 - R2) / -h);
                break;
        case 3: zwsp = ((y0 + h - R1) / h) * ((z0 - R2) / -h);
                break;
        default: std::cout << "Fehler: kein lokaler Knoten" << std::endl;
    }
    return zwsp;
}

}//namespace Icarus
