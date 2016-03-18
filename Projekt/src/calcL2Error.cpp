#include "include/assemblefem.hpp"

namespace Icarus
{

//TODO TOCHECK Superkonvergenz??
double assembleFem::calcL2Error(mathfunction realsol, FullVector<double> calcsol)
{
    double l2norm(0.0);
    std::vector<double> quadrature_X(27), quadrature_Y(27), quadrature_Z(27), basis3d(27);
    quadrature_X = get_quadrature_xpoints();
    quadrature_Y = get_quadrature_ypoints();
    quadrature_Z = get_quadrature_zpoints();
    //Hilfsvektor (DoF-Manager)
    std::vector<int> a{0, 1, _nx, _nx+1, _nx*_ny, _nx*_ny+1, _nx*_ny+_nx, _nx*_ny+_nx+1};

    //Schleife ueber alle Elemente
    for (int x(0); x < _nx-1; ++x)
        for (int y(0); y < _ny-1; ++y)
            for (int z(0); z < _nz-1; ++z)
            {
                int e(z*_ny*_nx+y*_nx+x);
                double e_x(getx(e));
                double e_y(gety(e));
                double e_z(getz(e));
                //Quadratur
                for(int q(0); q < 27; q++)
                {
                    double X(quadrature_X[q] + e_x);
                    double Y(quadrature_Y[q] + e_y);
                    double Z(quadrature_Z[q] + e_z);
                    //Auswertung von Gausspunkt q von der Loesung
                    double zwsp(0.0);
                    for(int A(0); A < 8; ++A)
                    {
                        basis3d = evaluated_Basis3d(A);
                        zwsp += basis3d[q] * calcsol[e + a[A]];
                    }
                    l2norm += (realsol.eval(X, Y, Z) - zwsp) * (realsol.eval(X, Y, Z) - zwsp) * _weight[q];
                }
            }
    return sqrt(l2norm);
}

}//namespace Icarus
