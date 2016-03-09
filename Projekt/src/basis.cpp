//TODO TODISCUSS:
// immer direkt return ...; ? also zwsp und x0,y0,z0 umgehen

#include "include/assemblefem.hpp"

namespace Icarus
{

double assembleFem::evaluate_Basis3d(int e, int A, double X, double Y, double Z)
{
    switch(A)
    {
        case 0: return ((getx(e) + _h - X) / _h) * ((gety(e) + _h - Y) / _h) * ((getz(e) + _h - Z) / _h);
        case 1: return ((getx(e) - X) / -_h) * ((gety(e) + _h - Y) / _h) * ((getz(e) + _h - Z) / _h);
        case 2: return ((getx(e) - X) / -_h) * ((gety(e) - Y) / -_h) * ((getz(e) + _h - Z) / _h);
        case 3: return ((getx(e) + _h - X) / _h) * ((gety(e) - Y) / -_h) * ((getz(e) + _h - Z) / _h);
        case 4: return ((getx(e) + _h - X) / _h) * ((gety(e) + _h - Y) / _h) * ((getz(e) - Z) / -_h);
        case 5: return ((getx(e) - X) / -_h) * ((gety(e) + _h - Y) / _h) * ((getz(e) - Z) / -_h);
        case 6: return ((getx(e) - X) / -_h) * ((gety(e) - Y) / -_h) * ((getz(e) - Z) / -_h);
        case 7: return ((getx(e) + _h - X) / _h) * ((gety(e) - Y) / -_h) * ((getz(e) - Z) / -_h);
    }
    assert(A > 0 && A < 8);
    //TODO Kompiler warnt natuerlich wegen evtl nicht erreichtem return. was machen?
}


std::vector<double> assembleFem::evaluate_gradient_Basis3d(int e, int A, double X, double Y, double Z)
{
    double x0(getx(e) - X);
    double y0(gety(e) - Y);
    double z0(getz(e) - Z);

    switch(A)
    {
        case 0: return {- (y0 + _h) * (z0 + _h) / (_h*_h*_h),
                        - (x0 + _h) * (z0 + _h) / (_h*_h*_h),
                        - (x0 + _h) * (y0 + _h) / (_h*_h*_h)};
        case 1: return {(y0 + _h) * (z0 + _h) / (_h*_h*_h),
                        x0 * (z0 + _h) / (_h*_h*_h),
                        x0 * (y0 + _h) / (_h*_h*_h)};
        case 2: return {- y0 * (z0 + _h) / (_h*_h*_h),
                        - x0 * (z0 + _h) / (_h*_h*_h),
                        - x0 * y0 / (_h*_h*_h)};
        case 3: return {y0 * (z0 + _h) / (_h*_h*_h),
                        (x0 + _h) * (z0 + _h) / (_h*_h*_h),
                        (x0 + _h) * y0 / (_h*_h*_h)};
        case 4: return {(y0 + _h) * z0 / (_h*_h*_h),
                        (x0 + _h) * z0 / (_h*_h*_h),
                        (x0 + _h) * (y0 + _h) / (_h*_h*_h)};
        case 5: return {- (y0 + _h) * z0 / (_h*_h*_h),
                        - x0 * z0 / (_h*_h*_h),
                        - x0 * (y0 + _h) / (_h*_h*_h)};
        case 6: return {y0 * z0 / (_h*_h*_h),
                        x0 * z0 / (_h*_h*_h),
                        x0 * y0 / (_h*_h*_h)};
        case 7: return {- y0 * z0 / (_h*_h*_h),
                        - (x0 + _h) * z0 / (_h*_h*_h),
                        - (x0 + _h) * y0 / (_h*_h*_h)};
    }
    assert(A > 0 && A < 8);
    //TODO Kompiler warnt natuerlich wegen evtl nicht erreichtem return. was machen?
}

double assembleFem::evaluate_Basis2d_1(int e, int A, double R1, double R2)
{
    double x0(getx(e));
    double y0(gety(e));
    double zwsp(0.0);

    switch (A)
    {
        case 0: zwsp = ((x0 + _h - R1) / _h) * ((y0 + _h - R2) / _h);
                break;
        case 1: zwsp = ((x0 - R1) / -_h) * ((y0 + _h - R2) / _h);
                break;
        case 2: zwsp = ((x0 - R1) / -_h) * ((y0 - R2) / -_h);
                break;
        case 3: zwsp = ((x0 + _h - R1) / _h) * ((y0 - R2) / -_h);
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
        case 0: zwsp = ((x0 + _h - R1) / _h) * ((z0 + _h - R2) / _h);
                break;
        case 1: zwsp = ((x0 - R1) / -_h) * ((z0 + _h - R2) / _h);
                break;
        case 2: zwsp = ((x0 - R1) / -_h) * ((z0 - R2) / -_h);
                break;
        case 3: zwsp = ((x0 + _h - R1) / _h) * ((z0 - R2) / -_h);
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
        case 0: zwsp = ((y0 + _h - R1) / _h) * ((z0 + _h - R2) / _h);
                break;
        case 1: zwsp = ((y0 - R1) / -_h) * ((z0 + _h - R2) / _h);
                break;
        case 2: zwsp = ((y0 - R1) / -_h) * ((z0 - R2) / -_h);
                break;
        case 3: zwsp = ((y0 + _h - R1) / _h) * ((z0 - R2) / -_h);
                break;
        default: std::cout << "Fehler: kein lokaler Knoten" << std::endl;
    }
    return zwsp;
}

}//namespace Icarus
