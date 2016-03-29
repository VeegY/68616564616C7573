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

class assembleFem
{
public:
    assembleFem(double sh, int sx, int sy, int sz):
        _h(sh), _nx(sx), _ny(sy), _nz(sz), z(_nx*_ny), y(_nx),
        _weight(get_weight()), _weight_2d(get_weight_2d()),
        _quadpoints_3d_x(get_quadrature_xpoints()), _quadpoints_3d_y(get_quadrature_ypoints()), _quadpoints_3d_z(get_quadrature_zpoints()),
        _quadpoints_2d_1(get_quadrature_points_2d_1()), _quadpoints_2d_2(get_quadrature_points_2d_2())
    { }

    void assemble(DistEllpackMatrix<double>& Matrix, SlicedVector<double>& rhs, std::vector<char>& disc_points,
        mathfunction f=mathfunction(0), mathfunction g=mathfunction(0), mathfunction h=mathfunction(0)); // rechte Seite, Dirichlet, Neumann

    double calcL2Error(mathfunction realsol, FullVector<double>& calcsol);

private:
    void assemblyMatrixRow(int rowlength);

    double assemblyRHSLoad(mathfunction f=mathfunction(0));
    double assemblyRHSNeumann(int Ebene, bool rightbacktop, mathfunction g=mathfunction(0));

    double getx(size_t index);
    double gety(size_t index);
    double getz(size_t index);

    int setup_A(int row);

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
    std::vector<double> get_quadrature_points_2d_1();
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
    std::vector<char>& _disc_points;
};

}//namespace Icarus

#endif//__ASSEMBLEFEM_HPP_
