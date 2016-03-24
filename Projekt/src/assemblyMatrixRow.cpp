// TODO TOCHECK: Ist das noch richtig?
// Werden die grad_Basen richtig eingesetzt?

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

    _value.clear();
    _column.resize(length);
    _value.resize(length);
    std::vector<bool> Belegt(length, false);

    //TODO TOCHECK fatal bei mehreren Objekten der Klasse, die unterschiedliche nx,ny und nz haben, oder?
    /*static*/ std::vector<int> a{0, 1, y, 1+y, z, 1+z, y+z, 1+y+z}; //Hilfsvektor um auf die Ecken (global gezaehlt) eines Element zu kommen

//    std::vector<std::vector<double>> grad_Basis1(3, std::vector<double>(27)); // 27 Mal so viel Speicher fuer 8 Mal weniger evalgradbasis aufrufen...
//    std::vector<std::vector<double>> grad_Basis2(3, std::vector<double>(27));
    /*static*/ std::vector<std::vector<std::vector<double>>> grad_Basis(8, std::vector<std::vector<double>>(3, std::vector<double>(27)));
//    static bool calced(false);
//    if (!calced)
        for (int B(0); B < 8; B++)
            grad_Basis[B] = evaluated_gradient_Basis3d(B);

    double zwsp(0.0);

    //Schleife ueber alle Elemente
    for (int i(0); i < n; i++)
    {
        //Berechne Grad_Basis3d fuer die Ecke A[i] vom Referenzelement
//        grad_Basis1 = evaluated_gradient_Basis3d(_A[i]);

        //Schleife ueber alle Ecken vom Referenzelement
        for (int B(0); B < 8; B++)
        {
            //Systematisch alle Ecken des Referenzelemets
//            grad_Basis2 = evaluated_gradient_Basis3d(B);

            //Quadratur
            zwsp = 0.0;
            for (int q(0); q < 27; q++)
            {
                // Quadratursumme
//                zwsp += (grad_Basis1[0][q]*grad_Basis2[0][q] + grad_Basis1[1][q]*grad_Basis2[1][q]
//                    + grad_Basis1[2][q]*grad_Basis2[2][q]) * _weight[q];
                zwsp += (grad_Basis[_A[i]][0][q]*grad_Basis[B][0][q] + grad_Basis[_A[i]][1][q]*grad_Basis[B][1][q]
                    + grad_Basis[_A[i]][2][q]*grad_Basis[B][2][q]) * _weight[q];
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
            assert(abbr);
        }//Schleife ueber alle Ecken von e
    }//Schleife ueber alle Elemente

//TODO TOCHECK sinnvoll?
    for (int i(0); i<length; ++i)
        if (_value[i] < 1.0e-9 && _value[i] > -1.0e-9)
            _value[i] = 0.0;
}//assemblyMatrixRow()

}//namespace Icarus
