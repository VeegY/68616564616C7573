#ifndef __ASSEMBLEFEM_HPP_
#define __ASSEMBLEFEM_HPP_

#include "distellpackmatrix.hpp"
#include "slicedvector.hpp"

#include "mathfunction.hpp"

#include <cmath>
#include <vector>

namespace Icarus
{

class assembleFem
{
public:
    assembleFem(double sh, int sx, int sy, int sz):
        h(sh), Nx(sx), Ny(sy), Nz(sz), z(Nx*Ny), y(Nx)
    {
        _weight = get_weight(); //TODO weight Vektor direkt bei Initialisierung setzen?
        _weight_2d = get_weight_2d();
        double a(sqrt(0.6));
        _ax = {-a, -a, -a, -a, -a, -a, -a, -a, -a, 0, 0, 0, 0, 0, 0, 0, 0, 0, a, a, a, a, a, a, a, a, a};  //x-Koordinaten der Gauss-Quadraturpunkte auf [-1,1]
        _ay = {-a, -a, -a, 0, 0, 0, a, a, a, -a, -a, -a, 0, 0, 0, a, a, a, -a, -a, -a, 0, 0, 0, a, a, a}; //y-Koordinaten der Gauss-Quadraturpunkte auf [-1,1]
        _az = {-a, 0, a, -a, 0, a, -a, 0, a, -a, 0, a, -a, 0, a, -a, 0, a, -a, 0, a, -a, 0, a, -a, 0, a};   //z-Koordinaten der Gauss-Quadraturpunkte auf [-1,1]
    }
    void assemble(DistEllpackMatrix<double>& Matrix, SlicedVector<double>& rhs);

//private:
    void assemblyMatrixRow(std::vector<int>& e, std::vector<int>& A, std::vector<int>& column, std::vector<double>& value);
    double assemblyRHSLoad(std::vector<int>& e, std::vector<int>& A, mathfunction f=mathfunction(0));
    double assemblyRHSNeumann(std::vector<int>& e, std::vector<int>& A, int Ebene, mathfunction g=mathfunction(0));
    double evaluate_Basis3d(int e, int A, double X, double Y, double Z);
    std::vector<double> evaluate_gradient_Basis3d(int e, int A, double X, double Y, double Z);
    double evaluate_Basis2d(int e, int A, int type, double R1, double R2);
    double getx(size_t index);
    double gety(size_t index);
    double getz(size_t index);
    std::vector<double> get_weight();
    std::vector<double> get_weight_2d();
    std::vector<double> get_quadrature_xpoints(int e);
    std::vector<double> get_quadrature_ypoints(int e);
    std::vector<double> get_quadrature_zpoints(int e);
    std::vector<double> get_quadrature_xpoints_2d(int e);
    std::vector<double> get_quadrature_ypoints_2d(int e);
    std::vector<double> get_quadrature_zpoints_2d(int e);

    double h;
    int Nx, Ny, Nz;
    int z, y;
    std::vector<double> _weight, _weight_2d;
    std::vector<double> _ax, _ay, _az;
};

}//namespace Icarus

#endif//__ASSEMBLEFEM_HPP_
