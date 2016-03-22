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
        h(sh), Nx(sx), Ny(sy), Nz(sz), z(Nx*Ny), y(Nx) {};
    void assemble(DistEllpackMatrix<double>& Matrix, SlicedVector<double>& rhs);

private:
    void assemblyMatrixRow(std::vector<int>& e, std::vector<int>& A, std::vector<int>& column, std::vector<double>& value);
    double assemblyRHSLoad(std::vector<int>& e, std::vector<int>& A, mathfunction f=mathfunction(0));
    double assemblyRHSNeumann(std::vector<int>& e, std::vector<int>& A, int Ebene, mathfunction g=mathfunction(0));
    double evaluate_Basis3d(int e, int A, double X, double Y, double Z);
    std::vector<double> evaluate_gradient_Basis3d(int e, int A, double X, double Y, double Z);
    double evaluate_Basis2d(int e, int A, int type, double R1, double R2);
    double getx(size_t index, double h, size_t nx, size_t ny);
    double gety(size_t index, double h, size_t nx, size_t ny);
    double getz(size_t index, double h, size_t nx, size_t ny);
    std::vector<double> get_weight(double c, double d);
    void transformation(std::vector<double>& ai, double h, std::vector<double>& trans);
    std::vector<double> get_quadrature_xpoints(int e, double h, std::vector<double>& ax, std::vector<double>& trans);
    std::vector<double> get_quadrature_ypoints(int e, double h, std::vector<double>& ay, std::vector<double>& trans);
    std::vector<double> get_quadrature_zpoints(int e, double h, std::vector<double>& az, std::vector<double>& trans);

    double h;
    int Nx, Ny, Nz;
    int z, y;
};

}//namespace Icarus

#endif//__ASSEMBLEFEM_HPP_