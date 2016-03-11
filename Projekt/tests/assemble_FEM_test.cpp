#include "../src/include/assemblefem.hpp"
#include "../src/include/mathfunction.hpp"

#include "../src/include/slicedvector.hpp"
#include "../src/include/bicgstabsolver.hpp"
#include "../src/include/distellpackmatrix.hpp"

#include "../src/include/vtkwriter.hpp"
//#include "../src/include/discretizer.hpp"

#include <iostream>
#include <cmath>

int main()
{
//*
    // u = 0
    // f = 0
    // d = 0
    double h(0.1);
    int nx(11), ny(11), nz(11);
    Icarus::mathfunction u(0, 0.0);
    Icarus::mathfunction f(0, 0.0);
    Icarus::mathfunction g(0, 0.0);

    Icarus::DistEllpackMatrix<double> matrix(nx*ny*nz);

    Icarus::SlicedVector<double> rhs(nx*ny*nz);

    Icarus::SlicedVector<double> sol(nx*ny*nz);
    sol.clear(); sol.set_local(0, 0.1);

    Icarus::assembleFem assembler(h, nx, ny, nz);
    assembler.assemble(matrix, rhs, f, g);

    Icarus::BiCgStabSolver<Icarus::DistEllpackMatrix<double>> solver(matrix, rhs);
    solver.solve(sol);

    Icarus::FullVector<double> fullsol(sol);
    Icarus::vtkWriter writer("f=0_dirichlet=0", "f=0_dirichlet=0", nx, ny, nz, h, 1);
    writer.addPointDataToTimestep(fullsol, 0, "Potential");

    double l2norm(0.0);
    for (int x(0); x < nx; ++x)
        for (int y(0); y < ny; ++y)
            for (int z(0); z < nz; ++z)
                l2norm += (u.eval(x*h, y*h, z*h)-fullsol[z*ny*nx+y*nx+x])*(u.eval(x*h, y*h, z*h)-fullsol[z*ny*nx+y*nx+x]);
    l2norm = sqrt(l2norm);
    LOG_INFO("l2norm = ", l2norm);
//*/
/*
    // u = 0
    // f = 0
    // d = 0
    // n = 0
    double h(0.1);
    int nx(11), ny(11), nz(11);
    Icarus::mathfunction u(0, 0.0);
    Icarus::mathfunction f(0, 0.0);
    Icarus::mathfunction g(0, 0.0);
    Icarus::mathfunction n(0, 0.0);

    Icarus::DistEllpackMatrix<double> matrix(nx*ny*nz);

    Icarus::SlicedVector<double> rhs(nx*ny*nz);

    Icarus::SlicedVector<double> sol(nx*ny*nz);
    sol.clear(); sol.set_local(0, 0.1);

    Icarus::assembleFem assembler(h, nx, ny, nz);
    assembler.assemble(matrix, rhs, f, g, n);

    Icarus::BiCgStabSolver<Icarus::DistEllpackMatrix<double>> solver(matrix, rhs);
    solver.solve(sol);

    Icarus::FullVector<double> fullsol(sol);
    Icarus::vtkWriter writer("f=0_dirichlet=0", "f=0_dirichlet=0", nx, ny, nz, h, 1);
    writer.addPointDataToTimestep(fullsol, 0, "Potential");

    double l2norm(0.0);
    for (int x(0); x < nx; ++x)
        for (int y(0); y < ny; ++y)
            for (int z(0); z < nz; ++z)
                l2norm += (u.eval(x*h, y*h, z*h)-fullsol[z*ny*nx+y*nx+x])*(u.eval(x*h, y*h, z*h)-fullsol[z*ny*nx+y*nx+x]);
    l2norm = sqrt(l2norm);
    LOG_INFO("l2norm = ", l2norm);
*/
/*
    // u = (0.5-x)(0.5-y)(0.5-z)
    // f = 0
    // d = Fallunterscheidung...
    double h(0.1);
    int nx(11), ny(11), nz(11);
    Icarus::mathfunction u(1);
    Icarus::mathfunction f(2);
    Icarus::mathfunction g(3);

    Icarus::DistEllpackMatrix<double> matrix(nx*ny*nz);

    Icarus::SlicedVector<double> rhs(nx*ny*nz);

    Icarus::SlicedVector<double> sol(nx*ny*nz);
    sol.clear(); sol.set_local(0, 0.1);

    Icarus::assembleFem assembler(h, nx, ny, nz);
    assembler.assemble(matrix, rhs, f, g);

    Icarus::BiCgStabSolver<Icarus::DistEllpackMatrix<double>> solver(matrix, rhs);
    solver.solve(sol);

    Icarus::FullVector<double> fullsol(sol);
    Icarus::vtkWriter writer("u=(0.5-x)...", "u=(0.5-x)...", nx, ny, nz, h, 1);
    writer.addPointDataToTimestep(fullsol, 0, "Potential");

    double l2norm(0.0);
    for (int x(0); x < nx; ++x)
        for (int y(0); y < ny; ++y)
            for (int z(0); z < nz; ++z)
                l2norm += (u.eval(x*h, y*h, z*h)-fullsol[z*ny*nx+y*nx+x])*(u.eval(x*h, y*h, z*h)-fullsol[z*ny*nx+y*nx+x]);
    l2norm = sqrt(l2norm);
    LOG_INFO("l2norm = ", l2norm);
*/
/*
    // u = (0.5-x)(0.5-y)(0.5-z)
    // f = 0
    // g = fallunterscheidung...
    // n = fallunterscheidung...
    double h(0.1);
    int nx(11), ny(11), nz(11);
    Icarus::mathfunction u(1);
    Icarus::mathfunction f(2);
    Icarus::mathfunction g(3);
    Icarus::mathfunction n(4);

    Icarus::DistEllpackMatrix<double> matrix(nx*ny*nz);

    Icarus::SlicedVector<double> rhs(nx*ny*nz);

    Icarus::SlicedVector<double> sol(nx*ny*nz);
    sol.clear(); sol.set_local(0, 0.1);

    Icarus::assembleFem assembler(h, nx, ny, nz);
    assembler.assemble(matrix, rhs, f, g, n);

    Icarus::BiCgStabSolver<Icarus::DistEllpackMatrix<double>> solver(matrix, rhs);
    solver.solve(sol);

    Icarus::FullVector<double> fullsol(sol);
    Icarus::vtkWriter writer("u=(0.5-x)...", "u=(0.5-x)...", nx, ny, nz, h, 1);
    writer.addPointDataToTimestep(fullsol, 0, "Potential");

    double l2norm(0.0);
    for (int x(0); x < nx; ++x)
        for (int y(0); y < ny; ++y)
            for (int z(0); z < nz; ++z)
                l2norm += (u.eval(x*h, y*h, z*h)-fullsol[z*ny*nx+y*nx+x])*(u.eval(x*h, y*h, z*h)-fullsol[z*ny*nx+y*nx+x]);
    l2norm = sqrt(l2norm);
    LOG_INFO("l2norm = ", l2norm);
*/
/*
    // u = x*(1.0-x)*y*(1.0-y)*z*(1.0-z)
    // f = 2.0*(x*(1.0-x)*y*(1.0-y)+x*(1.0-x)*z*(1.0-z)+y*(1.0-y)*z*(1.0-z))
    // d = 0
    double h(0.1);
    int nx(11), ny(11), nz(11);
    Icarus::mathfunction u(1);
    Icarus::mathfunction f(2);
    Icarus::mathfunction g(3);

    Icarus::DistEllpackMatrix<double> matrix(nx*ny*nz);

    Icarus::SlicedVector<double> rhs(nx*ny*nz);

    Icarus::SlicedVector<double> sol(nx*ny*nz);
    sol.clear(); sol.set_local(0, 0.1);

    Icarus::assembleFem assembler(h, nx, ny, nz);
    assembler.assemble(matrix, rhs, f, g);

    Icarus::BiCgStabSolver<Icarus::DistEllpackMatrix<double>> solver(matrix, rhs);
    solver.solve(sol);

    Icarus::FullVector<double> fullsol(sol);
    Icarus::vtkWriter writer("u=(0.5-x)...", "u=(0.5-x)...", nx, ny, nz, h, 1);
    writer.addPointDataToTimestep(fullsol, 0, "Potential");

    double l2norm(0.0);
    for (int x(0); x < nx; ++x)
        for (int y(0); y < ny; ++y)
            for (int z(0); z < nz; ++z)
                l2norm += (u.eval(x*h, y*h, z*h)-fullsol[z*ny*nx+y*nx+x])*(u.eval(x*h, y*h, z*h)-fullsol[z*ny*nx+y*nx+x]);
    l2norm = sqrt(l2norm);
    LOG_INFO("l2norm = ", l2norm);
*/

    return 0;
}
