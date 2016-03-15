#ifndef __ASSEMBLEFEM_HPP_
#define __ASSEMBLEFEM_HPP_

#include "distellpackmatrix.hpp"
#include "slicedvector.hpp"

#include "mathfunction.hpp"

#include <cmath>
#include <vector>
#include <cassert>

namespace Icarus
{

class assembleFem
{
public:
    assembleFem(double sh, int sx, int sy, int sz):
        _h(sh), _nx(sx), _ny(sy), _nz(sz), z(_nx*_ny), y(_nx)
    {
        _weight = get_weight(); //TODO weight Vektor direkt bei Initialisierung setzen?
        _weight_2d = get_weight_2d(); //TODO weight Vektor direkt bei Initialisierung setzen?
        double a(sqrt(0.6));
        _ax = {-a, -a, -a, -a, -a, -a, -a, -a, -a, 0, 0, 0, 0, 0, 0, 0, 0, 0, a, a, a, a, a, a, a, a, a}; //x-Koordinaten der Gauss-Quadraturpunkte auf [-1,1]
        _ay = {-a, -a, -a, 0, 0, 0, a, a, a, -a, -a, -a, 0, 0, 0, a, a, a, -a, -a, -a, 0, 0, 0, a, a, a}; //y-Koordinaten der Gauss-Quadraturpunkte auf [-1,1]
        _az = {-a, 0, a, -a, 0, a, -a, 0, a, -a, 0, a, -a, 0, a, -a, 0, a, -a, 0, a, -a, 0, a, -a, 0, a}; //z-Koordinaten der Gauss-Quadraturpunkte auf [-1,1]
    }

    void assemble(DistEllpackMatrix<double>& Matrix, SlicedVector<double>& rhs,
            mathfunction f=mathfunction(0), mathfunction g=mathfunction(0), mathfunction h=mathfunction(0)); // rechte Seite, Dirichlet, Neumann

private:
    void assemblyMatrixRow();

    double assemblyRHSLoad(mathfunction f=mathfunction(0));
    double assemblyRHSNeumann(int Ebene, int leftright, mathfunction g=mathfunction(0));

    double getx(size_t index);
    double gety(size_t index);
    double getz(size_t index);

    std::vector<double> evaluated_Basis3d(int A);
    std::vector<std::vector<double>> evaluated_gradient_Basis3d(int A);
    std::vector<double> get_weight();

    std::vector<double> get_quadrature_xpoints();
    std::vector<double> get_quadrature_ypoints();
    std::vector<double> get_quadrature_zpoints();

    std::vector<double> evaluated_Basis2d_1(int A);
    std::vector<double> evaluated_Basis2d_2(int A);
    std::vector<double> evaluated_Basis2d_3(int A);
    std::vector<double> get_weight_2d();
    std::vector<double> get_quadrature_xpoints_2d_1();
    std::vector<double> get_quadrature_ypoints_2d_1();
    std::vector<double> get_quadrature_xpoints_2d_2();
    std::vector<double> get_quadrature_zpoints_2d_2();
    std::vector<double> get_quadrature_ypoints_2d_3();
    std::vector<double> get_quadrature_zpoints_2d_3();

    double _h;
    int _nx, _ny, _nz;
    int z, y;
    std::vector<double> _weight, _weight_2d;
    std::vector<double> _ax, _ay, _az;
    std::vector<int> _e, _A;
    std::vector<int> _column;
    std::vector<double> _value;
};

}//namespace Icarus

#endif//__ASSEMBLEFEM_HPP_
