#include "include/assemblefem.hpp"

namespace Icarus
{

double assembleFem::assemblyRHSNeumann(std::vector<int>& e, std::vector<int>& A, int Ebene, mathfunction g)
{
    int n(e.size());
    double RHS(0.0);

    std::vector<double> R1(27), R2(27);

    //X-Y-Ebene
    if(Ebene == 1)
    {
        for(int i(0); i<n; i++)
        {
            //Quadraturpunkte fuer X-Y-Ebene
            R1 = get_quadrature_xpoints_2d_1(e[i]);
            R2 = get_quadrature_ypoints_2d_1(e[i]);
            int nqp = R1.size();

            double e_z(getz(e[i]));
            for(int q(0); q<nqp; q++)
                RHS += evaluate_Basis2d_1(e[i], A[i], R1[q], R2[q]) * g.eval(R1[q], R2[q], e_z) * _weight_2d[q];
        }
    }
    //X-Z-Ebene
    else if(Ebene == 2)
    {
        for(int i(0); i<n; i++)
        {
            //Quadraturpunkte fuer X-Z-Ebene
            R1 = get_quadrature_xpoints_2d_2(e[i]);
            R2 = get_quadrature_zpoints_2d_2(e[i]);
            int nqp(R1.size());

            double e_y(gety(e[i]));
            for(int q(0); q<nqp; q++)
                RHS += evaluate_Basis2d_2(e[i], A[i], R1[q], R2[q]) * g.eval(R1[q], e_y, R2[q]) * _weight_2d[q];
        }
    }
    //Y-Z-Ebene
    else if(Ebene == 3)
    {
        for(int i(0); i<n; i++)
        {
            //Quadraturpunkte fuer Y-Z-Ebene
            R1 = get_quadrature_ypoints_2d_3(e[i]);
            R2 = get_quadrature_zpoints_2d_3(e[i]);
            int nqp(R1.size());

            double e_x(getx(e[i]));
            for(int q(0); q<nqp; q++)
                RHS += evaluate_Basis2d_3(e[i], A[i], R1[q], R2[q]) * g.eval(e_x, R1[q], R2[q]) * _weight_2d[q];
        }
    }
    return RHS;
}//assemblyRHSNeumann

}//namespace Icarus


