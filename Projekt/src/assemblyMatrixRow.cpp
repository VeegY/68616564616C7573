//#include <vector>
//#include "quadratur.hpp"
//#include "basis.hpp"

#include "include/assemblefem.hpp"

namespace Icarus
{
//TODO TOCHECK e sollte size_t sein, anstatt von nur int
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
    //TODO TOCHECK Diese Vektoren sollte hier nur 0 Eintraege haben. Bin mir gerade nicht sicher, ob das der Fall ist
    column.clear();
    value.clear();
    column.resize(length);
    value.resize(length);

    std::vector<int> a(8); //Hilfsvektor um auf die Ecken (global gezaehlt) eines Element zu kommen
    a[0]=0; a[1]=1; a[2]= 1+y; a[3]=y; a[4]= z; a[5]=1+z; a[6]= 1+y+z; a[7]=y+z;

    std::vector<double> grad_Basis1(3);
    std::vector<double> grad_Basis2(3);

    double zwsp(0.0);

    //Fuer Quadratur
    //TODO TOCHECK changed 02-24-16
    std::vector<double> X(27), Y(27), Z(27);
    std::vector<double> transx(27);
    std::vector<double> transy(27);
    std::vector<double> transz(27);
    transformation(_ax, transx);
    transformation(_ay, transy);
    transformation(_az, transz);
    //TODO TOCHECK changed 02-24-16

    //Schleife ueber alle Elemente
    for(int i = 0; i < n; i++)
    {
        //getQuadrature(e[i], "Name") = X, Y, Z, weigth;
        //TODO TOCHECK changed 02-24-16
        X = get_quadrature_xpoints(e[i], transx);
        Y = get_quadrature_ypoints(e[i], transy);
        Z = get_quadrature_zpoints(e[i], transz);
        //TODO TOCHECK changed 02-24-16

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
            //Neuer Versuch.

            bool abbr(false);
            for(int j=0; (j<length) && (!abbr); j++)
            {
                if(column[j] == e[i] + a[B])
                {
std::cout << "i, B, j: " << i << ", " << B << ", " << j << " -> ";
std::cout << "set value[" << j << "] += " << zwsp << std::endl;
                    value[j] += zwsp;
                    abbr = true;
                }
                else if(j>0 && column[j] == 0)
                {
std::cout << "i, B, j: " << i << ", " << B << ", " << j << " -> ";
std::cout << "set value[" << j << "] += " << zwsp << ", column[" << j << "] = " << e[i]+a[B] << std::endl;
                    column[j] = e[i] + a[B];
                    value[j] += zwsp;
                    abbr = true;
                }
            }


            //Alte Variante. Hoechstwahrscheinlich Falsch
            /*
            int j=0;
            bool abbr = false;
            while(abbr==false && j<length)  // TODO TOCHECK und statt oder!?
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
            */
        }
    }

}

}//namespace Icarus
