//#include <vector>
//#include "quadratur.hpp"
//#include "basis.hpp"

#include "include/assemblefem.hpp"

namespace Icarus
{

void assembleFem::assemblyMatrixRow(std::vector<int>& e, std::vector<int>& A, std::vector<int>& column, std::vector<double>& value)
{
    //TODO sinnvoll?
    int z(Nx*Ny);
    int y(Nx);

    int n = e.size();
    int length(0);
    switch(n)
    {
    case 1:     length = 8;
                break;
    case 2:     length = 12;
                break;
    case 4:     length = 18;
                break;
    case 8:     length = 27;
                break;
    }

    column.resize(length);
    value.resize(length);

    std::vector<int> a(8);
    //std::vector<double> value(27);
    //std::vector<int> column(27);

    std::vector<double> grad_Basis1(3);
    std::vector<double> grad_Basis2(3);
    //std::vector<double> Zeile(Nx*Ny*Nz);

    double zwsp;

    a[0]=0; a[1]=1; a[2]= 1+y; a[3]=y; a[4]= z; a[5]=1+z; a[6]= 1+y+z; a[7]=y+z;

    //TODO TOCHECK changed 02-24-16
    std::vector<double> X(27), Y(27), Z(27), weight(27);
    double aa(sqrt(0.6));
    std::vector<double> ax{-aa, -aa, -aa, -aa, -aa, -aa, -aa, -aa, -aa, 0, 0, 0, 0, 0, 0, 0, 0, 0, aa, aa, aa, aa, aa, aa, aa, aa, aa};  //x-Koordinaten der Gauss-Quadraturpunkte auf [-1,1]
    std::vector<double> ay{-aa, -aa, -aa, 0.5, 0, 0, aa, aa, aa, -aa, -aa, -aa, 0, 0, 0, aa, aa, aa, -aa, -aa, -aa, 0, 0, 0.25, aa, aa, aa}; //y-Koordinaten der Gauss-Quadraturpunkte auf [-1,1]
    std::vector<double> az{-aa, 0, aa, -aa, 0, aa, -aa, 0, aa, -aa, 0, aa, -aa, 0, aa, -aa, 0, aa, -aa, 0, aa, -aa, 0, aa, -aa, 0, aa};   //z-Koordinaten der Gauss-Quadraturpunkte auf [-1,1]
    std::vector<double> trans(27);
    //TODO getweight fehlt noch
    //TODO TOCHECK changed 02-24-16

    //Schleife ueber alle Elemente
    for(int i = 0; i < n; i++)
    {
        //getQuadrature(e[i], "Name") = X, Y, Z, weigth;
        //TODO TOCHECK changed 02-24-16
        X = get_quadrature_xpoints(e[i], h, ax, trans);
        Y = get_quadrature_xpoints(e[i], h, ay, trans);
        Z = get_quadrature_xpoints(e[i], h, az, trans);
        //TODO TOCHECK changed 02-24-16

        int nqp = X.size();
        //Schleife ueber alle Ecken von e
        for(int B = 0; B < 8; B++)
        {
            //Quadratur
            zwsp = 0;
            for(int q = 0; q<nqp; q++)
            {
                grad_Basis1 = evaluate_gradient_Basis3d(e[i], A[i], X[q], Y[q], Z[q]);
                grad_Basis2 = evaluate_gradient_Basis3d(e[i], B, X[q], Y[q], Z[q]);
                //zwsp += grad_Basis1.dot(grad_Basis2) * weight[q];
                zwsp += (grad_Basis1[0]*grad_Basis2[0] + grad_Basis1[1]*grad_Basis2[1]
                    + grad_Basis1[2]*grad_Basis2[2]) * weight[q];
            }

            int j=0;
            bool abbr = false;
            while(abbr==false||j<length)
            {
                if(e[i] + a[B] == 0)
                {
                    column[j] = e[i] + a[B];
                    value[j] += zwsp;
                    abbr = true;
                }
                else if(e[i] + a[B] == column[j])
                {
                    value[j] += zwsp;
                    abbr = true;
                }
                j++;
            }
            //Zeile[e[i] + a[B]] += zwsp;
        }
    }

    //Zeile zu value und column

}

}//namespace Icarus