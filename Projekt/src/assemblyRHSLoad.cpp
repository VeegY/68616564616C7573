#include "include/assemblefem.hpp"

namespace Icarus
{

double assembleFem::assemblyRHSLoad(mathfunction f)
{
    int n(_e.size());
    double RHS(0.0);

//    std::vector<double> X(27), Y(27), Z(27), Basis3d(27);
    std::vector<double> Basis3d(27);

    for (int i(0); i < n; i++)
    {
        double e_x(getx(_e[i]));
        double e_y(gety(_e[i]));
        double e_z(getz(_e[i]));

        Basis3d = evaluated_Basis3d(_A[i]);

        //Zum auswerten von f, translatiere die Gauspunkte zum Element e[i]
        for (int q(0); q < 27; q++)
            RHS += Basis3d[q] * f.eval(_quadpoints_3d_x[q] + e_x, _quadpoints_3d_y[q] + e_y, _quadpoints_3d_z[q] + e_z) * _weight[q];
    }
    return RHS;
}

}//namespace Icarus
