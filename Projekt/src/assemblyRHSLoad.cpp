#include "include/assemblefem.hpp"

namespace Icarus
{

double assembleFem::assemblyRHSLoad(mathfunction f)
{
    int n(_e.size());
    double RHS(0.0);

    std::vector<double> X(27), Y(27), Z(27), Basis3d(27);

    for (int i(0); i < n; i++)
    {
        X = get_quadrature_xpoints();
        Y = get_quadrature_ypoints();
        Z = get_quadrature_zpoints();
        double e_x(getx(_e[i]));
        double e_y(gety(_e[i]));
        double e_z(getz(_e[i]));

        Basis3d = evaluated_Basis3d(_A[i]);

        for (int q(0); q < 27; q++)
            //Zum auswerten von f, translatiere die Gauspunkte zum Element e[i]
            RHS += Basis3d[q] * f.eval(X[q] + e_x, Y[q] + e_y, Z[q] + e_z) * _weight[q];
    }
    return RHS;
}

}//namespace Icarus
