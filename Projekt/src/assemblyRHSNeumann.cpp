#include "include/assemblefem.hpp"

namespace Icarus
{

double assembleFem::assemblyRHSNeumann(int Ebene, bool rightbacktop, mathfunction g)
{
    assert(Ebene == 1 || Ebene == 2 || Ebene == 3);

    int n(_e.size());
    double RHS(0.0);

    std::vector<double> Basis2d(9);

    //X-Y-Ebene
    if(Ebene == 1)
    {
        for(int i(0); i<n; i++)
        {
            double e_x(getx(_e[i]));
            double e_y(gety(_e[i]));
            double e_z((rightbacktop ? getz(_e[i]) + _h : getz(_e[i])));
            Basis2d = (rightbacktop ? evaluated_Basis2d_1(_A[i] - 4) : evaluated_Basis2d_1(_A[i]));

            for(int q(0); q < 9; q++)
            {
                //Zum auswerten von g, translatiere die Gauspunkte zum Element e[i]
                if (rightbacktop)
                    RHS += Basis2d[q] * g.eval(_quadpoints_2d_1[q] + e_x, _quadpoints_2d_2[q] + e_y, e_z, Ebene) * _weight_2d[q];
                else
                    RHS -= Basis2d[q] * g.eval(_quadpoints_2d_1[q] + e_x, _quadpoints_2d_2[q] + e_y, e_z, Ebene) * _weight_2d[q];
            }
        }
    }
    //X-Z-Ebene
    else if(Ebene == 2)
    {
        for(int i(0); i<n; i++)
        {
            double e_x(getx(_e[i]));
            double e_y((rightbacktop ? gety(_e[i]) + _h : gety(_e[i])));
            double e_z(getz(_e[i]));
            Basis2d = (rightbacktop ? evaluated_Basis2d_2(_A[i] - 2) : evaluated_Basis2d_2(_A[i]));

            for(int q(0); q < 9; q++)
            {
                //Zum auswerten von g, translatiere die Gauspunkte zum Element e[i]
                if (rightbacktop)
                    RHS += Basis2d[q] * g.eval(_quadpoints_2d_1[q] + e_x, e_y, _quadpoints_2d_2[q] + e_z, Ebene) * _weight_2d[q];
                else
                    RHS -= Basis2d[q] * g.eval(_quadpoints_2d_1[q] + e_x, e_y, _quadpoints_2d_2[q] + e_z, Ebene) * _weight_2d[q];
            }
        }
    }
    //Y-Z-Ebene
    else if(Ebene == 3)
    {
        for(int i(0); i<n; i++)
        {
            double e_x((rightbacktop ? getx(_e[i]) + _h : getx(_e[i])));
            double e_y(gety(_e[i]));
            double e_z(getz(_e[i]));
            Basis2d = (rightbacktop ? evaluated_Basis2d_3(_A[i] - 1) : evaluated_Basis2d_3(_A[i]));

            for(int q(0); q < 9; q++)
            {
                //Zum auswerten von g, translatiere die Gauspunkte zum Element e[i]
                if (rightbacktop)
                    RHS += Basis2d[q] * g.eval(e_x, _quadpoints_2d_1[q] + e_y, _quadpoints_2d_2[q] + e_z, Ebene) * _weight_2d[q];
                else
                    RHS -= Basis2d[q] * g.eval(e_x, _quadpoints_2d_1[q] + e_y, _quadpoints_2d_2[q] + e_z, Ebene) * _weight_2d[q];
            }
        }
    }

    return RHS;
}//assemblyRHSNeumann

}//namespace Icarus
