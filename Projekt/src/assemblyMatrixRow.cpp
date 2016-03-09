#include "include/assemblefem.hpp"

namespace Icarus
{
void assembleFem::assemblyMatrixRow()
{
    int n(_e.size()); //Anzahl zu betrachtender Elemente
    int length(0);
    switch(n) //Anpassen der Laenge der Ausgabe
    {
    case 1: length = 8;
            break;
    case 2: length = 12;
            break;
    case 4: length = 18;
            break;
    case 8: length = 27;
            break;
    default: assert(_e.size() == 1 || _e.size() == 2 || _e.size() == 4 || _e.size() == 8);
    }

    //column.clear();
    _value.clear();
    _column.resize(length);
    _value.resize(length);
    std::vector<bool> Belegt(length, false);

    //std::vector<int> a(8); //Hilfsvektor um auf die Ecken (global gezaehlt) eines Element zu kommen
    //a[0]=0; a[1]=1; a[2]= 1+y; a[3]=y; a[4]= z; a[5]=1+z; a[6]= 1+y+z; a[7]=y+z;
    static std::vector<int> a{0, 1, 1+y, y, z, 1+z, 1+y+z, y+z};

    std::vector<std::vector<double>> grad_Basis1(27, std::vector<double>(3)); // 27 Mal so viel Speicher fuer 8 Mal weniger evalgradbasis aufrufen...
    std::vector<double> grad_Basis2(3);

    double zwsp(0.0);

    //Fuer Quadratur
    std::vector<double> X(27), Y(27), Z(27);

    //Schleife ueber alle Elemente
    for (int i(0); i < n; i++)
    {
        //Quadraturpunkte festlegen
        X = get_quadrature_xpoints(_e[i]);
        Y = get_quadrature_ypoints(_e[i]);
        Z = get_quadrature_zpoints(_e[i]);
        for (int q(0); q < 27; q++)
            grad_Basis1[q] = evaluate_gradient_Basis3d(_e[i], _A[i], X[q], Y[q], Z[q]);

        //Schleife ueber alle Ecken von e
        for (int B(0); B < 8; B++)
        {
            //Quadratur
            zwsp = 0.0;
            for (int q(0); q < 27; q++)
            {
                //Berechne Grad_Basis3d fuer die Ecke A[i] von Element e[i]. Dies sollte immer der momentane Raumpunkt ('Zeile') sein.
                //grad_Basis1 = evaluate_gradient_Basis3d(e[i], A[i], X[q], Y[q], Z[q]);
                //Systematisch alle Ecken des Elemets e[i]
                grad_Basis2 = evaluate_gradient_Basis3d(_e[i], B, X[q], Y[q], Z[q]);

                //zwsp += grad_Basis1.dot(grad_Basis2) * weight[q]; (Also die Quadratursumme)
                zwsp += (grad_Basis1[q][0]*grad_Basis2[0] + grad_Basis1[q][1]*grad_Basis2[1]
                    + grad_Basis1[q][2]*grad_Basis2[2]) * _weight[q];
            }

            //Berechneter Wert an die richtige Stelle von Column und Value aufaddieren. Ich schätze, dass hier irgenwo der Fehler liegt.
            bool abbr(false);
            for(int j(0); (j<length) && (!abbr); j++)
            {
                if(!Belegt[j])
                {
                    _value[j] += zwsp;
                    _column[j] = _e[i] + a[B];
                    abbr = true;
                    Belegt[j] = true;
                }
                else if(_column[j] == _e[i] + a[B])
                {
                    _value[j] += zwsp;
                    abbr = true;
                }
            }
        }//Schleife ueber alle Ecken von e
    }//Schleife ueber alle Elemente
}//assemblyMatrixRow()

}//namespace Icarus
