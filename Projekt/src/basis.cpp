//#include "include/basis.hpp"
#include "include/assemblefem.hpp"

namespace Icarus
{

double assembleFem::evaluate_Basis3d(int e, int A, double X, double Y, double Z)
{
    //double x0(getXcoordinate(e));
    //double y0(getYcoordinate(e));
    //double z0(getZcoordinate(e));
    double x0(getx(e));   //Setzen von x0, y0, z0 und h nur zum Testen, ersetze durch auskommentiertes getcoordinate
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
    //double x0(getXcoordinate(e));
    //double y0(getYcoordinate(e));
    //double z0(getZcoordinate(e));
    double x0(getx(e));   //Setzen von x0, y0, z0 und h nur zum Testen, ersetze durch auskommentiertes getcoordinate
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

double assembleFem::evaluate_Basis2d(int e, int A, int type, double R1, double R2)
{
    //double x0(getXcoordinate(e));
    //double y0(getYcoordinate(e));
    //double z0(getZcoordinate(e));
    double x0(getx(e));   //Setzen von x0, y0, z0 und h nur zum Testen, ersetze durch auskommentiertes getcoordinate
    double y0(gety(e));
    double z0(getz(e));
    double zwsp(0.0);

    switch(type)
    {
        case 0: if(A==0)
                {
                    zwsp = ((x0 + h - R1) / h) * ((y0 + h - R2) / h);
                }
                else if(A==1)
                {
                    zwsp = ((x0 - R1) / -h) * ((y0 + h - R2) / h);
                }
                else if(A==2)
                {
                    zwsp = ((x0 - R1) / -h) * ((y0 - R2) / -h);
                }
                else if(A==3)
                {
                    zwsp = ((x0 + h - R1) / h) * ((y0 - R2) / -h);
                }
                break;
        case 1: if(A==0)
                {
                    zwsp = ((x0 + h - R1) / h) * ((z0 + h - R2) / h);
                }
                else if(A==1)
                {
                    zwsp = ((x0 - R1) / -h) * ((z0 + h - R2) / h);
                }
                else if(A==2)
                {
                    zwsp = ((x0 - R1) / -h) * ((z0 - R2) / -h);
                }
                else if(A==3)
                {
                    zwsp = ((x0 + h - R1) / h) * ((z0 - R2) / -h);
                }
                break;
        case 2: if(A==0)
                {
                    zwsp = ((y0 + h - R1) / h) * ((z0 + h - R2) / h);
                }
                else if(A==1)
                {
                    zwsp = ((y0 - R1) / -h) * ((z0 + h - R2) / h);
                }
                else if(A==2)
                {
                    zwsp = ((y0 - R1) / -h) * ((z0 - R2) / -h);
                }
                else if(A==3)
                {
                    zwsp = ((y0 + h - R1) / h) * ((z0 - R2) / -h);
                }
                break;
        default: std::cout << "Fehler: keine gültige Ebene" << std::endl;
    }
    return zwsp;
}

}//namespace Icarus
