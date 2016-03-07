#include "include/assemblefem.hpp"

namespace Icarus
{
void assembleFem::assemblyMatrixRow(std::vector<int>& e, std::vector<int>& A, std::vector<int>& column, std::vector<double>& value)
{
    int n = e.size(); //Anzahl zu betrachtender Elemente
    int length(0);
    switch(n) //Anpassen der Laenge der Ausgabe
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

    column.clear();
    value.clear();
    column.resize(length);
    value.resize(length);
    std::vector<bool> Belegt(length, false);

    std::vector<int> a(8); //Hilfsvektor um auf die Ecken (global gezaehlt) eines Element zu kommen
    a[0]=0; a[1]=1; a[2]= 1+y; a[3]=y; a[4]= z; a[5]=1+z; a[6]= 1+y+z; a[7]=y+z;

    std::vector<double> grad_Basis1(3);
    std::vector<double> grad_Basis2(3);

    double zwsp(0.0);

    //Fuer Quadratur
    std::vector<double> X(27), Y(27), Z(27);

    //Schleife ueber alle Elemente
    for(int i = 0; i < n; i++)
    {
        //Quadraturpunkte festlegen
        X = get_quadrature_xpoints(e[i]);
        Y = get_quadrature_ypoints(e[i]);
        Z = get_quadrature_zpoints(e[i]);
        int nqp = X.size();

        //Schleife ueber alle Ecken von e
        for(int B = 0; B < 8; B++)
        {
            //Quadratur
            zwsp = 0.0;
            for(int q = 0; q<nqp; q++)
            {
                //Berechne Grad_Basis3d fuer die Ecke A[i] von Element e[i]. Dies sollte immer der momentane Raumpunkt ('Zeile') sein.
                grad_Basis1 = evaluate_gradient_Basis3d(e[i], A[i], X[q], Y[q], Z[q]);
                //Systematisch alle Ecken des Elemets e[i]
                grad_Basis2 = evaluate_gradient_Basis3d(e[i], B, X[q], Y[q], Z[q]);

                //zwsp += grad_Basis1.dot(grad_Basis2) * weight[q]; (Also die Quadratursumme)
                zwsp += (grad_Basis1[0]*grad_Basis2[0] + grad_Basis1[1]*grad_Basis2[1]
                    + grad_Basis1[2]*grad_Basis2[2]) * _weight[q];
            }

            //Berechneter Wert an die richtige Stelle von Column und Value aufaddieren. Ich schätze, dass hier irgenwo der Fehler liegt.
            bool abbr(false);
            for(int j=0; (j<length) && (!abbr); j++)
            {
                if(!Belegt[j])
                {
//std::cout << "i, B, j: " << i << ", " << B << ", " << j << " -> ";
//std::cout << "set value[" << j << "] += " << zwsp << ", column[" << j << "] = " << e[i]+a[B] << std::endl;
                    value[j] += zwsp;
                    column[j] = e[i] + a[B];
                    abbr = true;
                    Belegt[j] = true;
                }
                else if((column[j] == e[i] + a[B]))
                {
//std::cout << "i, B, j: " << i << ", " << B << ", " << j << " -> ";
//std::cout << "set value[" << j << "] += " << zwsp << std::endl;
                    value[j] += zwsp;
                    abbr = true;
                }
            }
        }
    }
}

}//namespace Icarus

