#ifndef __ASSEMBLEFEM_HPP_
#define __ASSEMBLEFEM_HPP_

#include "distellpackmatrix.hpp"
#include "slicedvector.hpp"
#include "fullvector.hpp"

#include "mathfunction.hpp"

#include <cmath>
#include <vector>
#include <cassert>

namespace Icarus
{
/**
*   \brief Assemblierung der Matrix und der rechten Seite für ein LGS mit einer finiten Elemente Methode.
*   In der Klasse enthalten ist eine 3D und eine 2D Quadratur, sowie die 3D/2D Basisfunktionen.
*   Mit den assembly-Funktionen werden die LGS-Einträge Zeilenweise berechnet, sodass die Assemblierung parallel laufen kann.
*   Um die Konvergenz zu überprüfen, ist zusätzlich eine Fehleberechnung bzgl. der L2-Norm vorhanden.
**/
class assembleFem
{
public:

    /**
    *   \brief Standardkonstruktor.
    *
    *   Da einige Werte (zum Beispiel die Auswertungen der Quadraturpunkte in den Basisfunktionen) mehrmals genutzt werden,
    *   werden diese global abgespeichert, sodass die Berechnungszeit verringert wird.
    *
    *   \param sh  Länge einer Schrittweite
    *   \param sx  Anzahl Stützstellen in Richtung der X-Achse
    *   \param sy  Anzahl Stützstellen in Richtung der Y-Achse
    *   \param sz  Anzahl Stützstellen in Richtung der Z-Achse
    *
    **/


    assembleFem(double sh, int sx, int sy, int sz):
        _h(sh), _nx(sx), _ny(sy), _nz(sz), z(_nx*_ny), y(_nx),
        _weight(get_weight()), _weight_2d(get_weight_2d()),
        _quadpoints_3d_x(get_quadrature_xpoints()), _quadpoints_3d_y(get_quadrature_ypoints()), _quadpoints_3d_z(get_quadrature_zpoints()),
        _quadpoints_2d_1(get_quadrature_points_2d_1()), _quadpoints_2d_2(get_quadrature_points_2d_2())
    { }

    /**
    *   \brief LGS-Assemblierung.
    *   Speichert die berechneten Werte in der vorher angelegten Ellpack-Matrix und den Sliced-Vector.
    *
    *   \param Matrix      Matrix des LGS
    *   \param rhs         Rechte Seite des LGS
    *   \param disc_points Typ der einzelnen Raumpunkte (aus dem discretizer)
    *   \param f           rechte Seite vom Problem -laplcae(u) = f
    *   \param g           Fuktion für die Dirichlet-Randbedingungen
    *   \param h           Funktion für die Neumann-Randbedingungen
    *
    **/
    void assemble(DistEllpackMatrix<double>& Matrix, SlicedVector<double>& rhs, std::vector<char>& disc_points,
        mathfunction f=mathfunction(0), mathfunction g=mathfunction(0), mathfunction h=mathfunction(0)); // rechte Seite, Dirichlet, Neumann

    /**
    *   \brief Fehlerberechnung bzgl. der L2-Norm
    *   falls die analytische Lösung des Problems bekannt ist,
    *   kann man mit dieser Funktion den Fehler u-u_h bzgl. der L2-Norm berechnen
    *
    *   \param realsol     analytische Lösung des Problems -laplace(u) = f
    *   \param calcsol     Lösung des LGS
    **/
    double calcL2Error(mathfunction realsol, FullVector<double>& calcsol);

private:

    /**
    *   \brief assembliert eine ganze Matrixzeile
    *
    *   \param rowlength   Anzahl der Einträge einer Zeile
    **/
    void assemblyMatrixRow(int rowlength);

    /**
    *   \brief berechnet das Integral der Lastfunktion. Dieser Wert wird zur assemblierung der rechten Seite genutzt.
    *
    *   \param f   rechte Seite von -laplace(u) = f
    **/
    double assemblyRHSLoad(mathfunction f=mathfunction(0));

    /**
    *   \brief Berechnet das Integral der Neumannfunktion. Dieser Wert wird zur assemblierung der rechten Seite genutzt.
    *
    *   \param Ebene           Bestimmt den 'Typ' der Ebene. 1->XY-Ebene, 2->XZ-Ebene, 3->YZ-Ebene
    *   \param rightbacktop    Bestimmt die Richtung, in der der Gradient zeigt
    *   \param g               Die Funktion der Neumann-Randbedingung
    **/
    double assemblyRHSNeumann(int Ebene, bool rightbacktop, mathfunction g=mathfunction(0));

    /**
    *   \brief Berechnet die X-Koordinate eines Raumpunktes
    *   \param index   Der globale Index des Raumpunktes
    **/
    double getx(size_t index);

    /**
    *   \brief Berechnet die Y-Koordinate eines Raumpunktes
    *   \param index   Der globale Index des Raumpunktes
    **/
    double gety(size_t index);

    /**
    *   \brief Berechnet die Z-Koordinate eines Raumpunktes
    *   \param index   Der globale Index des Raumpunktes
    **/
    double getz(size_t index);

    /**
    *   \brief Konfiguriert _A
    *   \param row         globaler Index vom Raumpunkt
    *   \param disc_points Typus der Raumpunkte
    **/
    int setup_A(int row, std::vector<char>& disc_points);

    /**
    *   \brief Konfiguriert _e
    *   \param row         globaler Index vom Raumpunkt
    **/
    void setup_e(int row);

    //TODO
    //void setup_neumann(int row, int plane, bool rightbacktop, std::vector<char>& disc_points);

    /**
    *   \brief Man wählt eine der 3D-Basen aus und erhält die ausgewerteten Quadraturpunkte
    *
    *   \param A   Bestimmt die Basis
    *   \return Auswertung in den Quadraturpunkten
    **/
    std::vector<double> evaluated_Basis3d(int A);

    /**
    *   \brief Man wählt einen der Gradienten der 3D-Basen aus und erhält die ausgewerteten Quadraturpunkte
    *
    *   \param A   Bestimmt die Basis
    *   \return Auswertung in den Quadraturpunkten
    **/
    std::vector<std::vector<double>> evaluated_gradient_Basis3d(int A);

    /**
    *   \brief Gewichte der 3D-Quadratur
    **/
    std::vector<double> get_weight();

    /**
    *   \brief X Koordinaten der 3D-Quadratur
    **/
    std::vector<double> get_quadrature_xpoints();

    /**
    *   \brief Y Koordinaten der 3D-Quadratur
    **/
    std::vector<double> get_quadrature_ypoints();

    /**
    *   \brief Z Koordinaten der 3D-Quadratur
    **/
    std::vector<double> get_quadrature_zpoints();

    /**
    *   \brief Man wählt einen der 2D-Basen aus und erhält die ausgewerteten Quadraturpunkte. Dies ist für die XY-Ebene
    **/
    std::vector<double> evaluated_Basis2d_1(int A);

    /**
    *   \brief Man wählt einen der 2D-Basen aus und erhält die ausgewerteten Quadraturpunkte. Dies ist für die XZ-Ebene
    **/
    std::vector<double> evaluated_Basis2d_2(int A); //brauchen wir überhaupt noch die 3 verschiedenen Fälle?

    /**
    *   \brief Man wählt einen der 2D-Basen aus und erhält die ausgewerteten Quadraturpunkte. Dies ist für die YZ-Ebene
    **/
    std::vector<double> evaluated_Basis2d_3(int A);

    /**
    *   \brief Gewichte der 2D-Quadratur
    **/
    std::vector<double> get_weight_2d();

    /**
    *   \brief Koordinaten der 2D-Quadratur in Richung 1 (R1)
    **/
    std::vector<double> get_quadrature_points_2d_1();

    /**
    *   \brief Koordinaten der 2D-Quadratur in Richung 2 (R2)
    **/
    std::vector<double> get_quadrature_points_2d_2();

    double _h;
    int _nx, _ny, _nz;
    int z, y;
    std::vector<double> _weight, _weight_2d;
    std::vector<double> _quadpoints_3d_x, _quadpoints_3d_y, _quadpoints_3d_z;
    std::vector<double> _quadpoints_2d_1, _quadpoints_2d_2;
    std::vector<int> _e, _A;
    std::vector<int> _column;
    std::vector<double> _value;
};

}//namespace Icarus

#endif//__ASSEMBLEFEM_HPP_
