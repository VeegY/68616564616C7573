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
//    double h(1.0/2.0);
    double h(1.0/40.0);
    int nx(41), ny(41), nz(41);
//    int nx(3), ny(3), nz(3);

    // ***** ***** ***** ***** TEST 0 DIRICHLET ***** ***** ***** ***** //
    // u = 0
    // f = 0
    // d = 0
    Icarus::mathfunction u_0(0, 0.0);
    Icarus::mathfunction f_0(0, 0.0);
    Icarus::mathfunction d_0(0, 0.0);
/*
    Icarus::DistEllpackMatrix<double> matrix_0d(nx*ny*nz);
    Icarus::SlicedVector<double> rhs_0d(nx*ny*nz);
    Icarus::SlicedVector<double> sol_0d(nx*ny*nz);
    sol_0d.clear(); sol_0d.set_local(0, 0.1);

    Icarus::assembleFem assembler_0d(h, nx, ny, nz);
    assembler_0d.assemble(matrix_0d, rhs_0d, f_0, d_0);

    Icarus::BiCgStabSolver<Icarus::DistEllpackMatrix<double>> solver_0d(matrix_0d, rhs_0d);
    solver_0d.solve(sol_0d);

    Icarus::FullVector<double> fullsol_0d(sol_0d);
    Icarus::vtkWriter writer_0d("f=0_dirichlet=0", "f=0_dirichlet=0", nx, ny, nz, h, 1);
    writer_0d.addPointDataToTimestep(fullsol_0d, 0, "Potential");

    double l2norm_0d(0.0);
    for (int x(0); x < nx; ++x)
        for (int y(0); y < ny; ++y)
            for (int z(0); z < nz; ++z)
                l2norm_0d += (u_0.eval(x*h, y*h, z*h)-fullsol_0d[z*ny*nx+y*nx+x])*(u_0.eval(x*h, y*h, z*h)-fullsol_0d[z*ny*nx+y*nx+x]);
    l2norm_0d = sqrt(l2norm_0d);
    LOG_INFO("l2norm_0d = ", l2norm_0d);
*/
///*
    // ***** ***** ***** ***** TEST 0 NEUMANN ***** ***** ***** ***** //
    // u = 0
    // f = 0
    // d = 0
    // n = 0
    Icarus::mathfunction n_0(0, 0.0);

    Icarus::DistEllpackMatrix<double> matrix_0n(nx*ny*nz);
    Icarus::SlicedVector<double> rhs_0n(nx*ny*nz);
    Icarus::SlicedVector<double> sol_0n(nx*ny*nz);
    sol_0n.clear(); sol_0n.set_local(0, 0.1);

    Icarus::assembleFem assembler_0n(h, nx, ny, nz);
    assembler_0n.assemble(matrix_0n, rhs_0n, f_0, d_0, n_0);

    Icarus::BiCgStabSolver<Icarus::DistEllpackMatrix<double>> solver_0n(matrix_0n, rhs_0n);
    solver_0n.solve(sol_0n);

    Icarus::FullVector<double> fullsol_0n(sol_0n);
    Icarus::vtkWriter writer_0n("f=0_dirichlet=0", "f=0_dirichlet=0", nx, ny, nz, h, 1);
    writer_0n.addPointDataToTimestep(fullsol_0n, 0, "Potential");

    double l2norm_0n(0.0);
    for (int x(0); x < nx; ++x)
        for (int y(0); y < ny; ++y)
            for (int z(0); z < nz; ++z)
                l2norm_0n += (u_0.eval(x*h, y*h, z*h)-fullsol_0n[z*ny*nx+y*nx+x])*(u_0.eval(x*h, y*h, z*h)-fullsol_0n[z*ny*nx+y*nx+x]);
    l2norm_0n = sqrt(l2norm_0n);
    LOG_INFO("l2norm_0n = ", l2norm_0n);

//matrix_0n.print_local_data(std::cout);
//*/
    // ***** ***** ***** ***** TEST 1 DIRICHLET ***** ***** ***** ***** //
    // u = (0.5-x)(0.5-y)(0.5-z)
    // f = 0
    // d = Fallunterscheidung...
    Icarus::mathfunction u_1(1);
    Icarus::mathfunction f_1(2);
    Icarus::mathfunction d_1(3);
/*
    Icarus::DistEllpackMatrix<double> matrix_1d(nx*ny*nz);
    Icarus::SlicedVector<double> rhs_1d(nx*ny*nz);
    Icarus::SlicedVector<double> sol_1d(nx*ny*nz);
    sol_1d.clear(); sol_1d.set_local(0, 0.1);

    Icarus::assembleFem assembler_1d(h, nx, ny, nz);
    assembler_1d.assemble(matrix_1d, rhs_1d, f_1, d_1);

    Icarus::BiCgStabSolver<Icarus::DistEllpackMatrix<double>> solver_1d(matrix_1d, rhs_1d);
    solver_1d.solve(sol_1d);

    Icarus::FullVector<double> fullsol_1d(sol_1d);
    Icarus::vtkWriter writer_1d("u=(0.5-x)...", "u=(0.5-x)...", nx, ny, nz, h, 1);
    writer_1d.addPointDataToTimestep(fullsol_1d, 0, "Potential");

    double l2norm_1d(0.0);
    for (int x(0); x < nx; ++x)
        for (int y(0); y < ny; ++y)
            for (int z(0); z < nz; ++z)
                l2norm_1d += (u_1.eval(x*h, y*h, z*h)-fullsol_1d[z*ny*nx+y*nx+x])*(u_1.eval(x*h, y*h, z*h)-fullsol_1d[z*ny*nx+y*nx+x]);
    l2norm_1d = sqrt(l2norm_1d);
    LOG_INFO("l2norm_1d = ", l2norm_1d);
*/
    // ***** ***** ***** ***** TEST 1 NEUMANN ***** ***** ***** ***** //
    // u = (0.5-x)(0.5-y)(0.5-z)
    // f = 0
    // g = Fallunterscheidung...
    // n = Fallunterscheidung...
    Icarus::mathfunction n_1(4);
/*
    Icarus::DistEllpackMatrix<double> matrix_1n(nx*ny*nz);
    Icarus::SlicedVector<double> rhs_1n(nx*ny*nz);
    Icarus::SlicedVector<double> sol_1n(nx*ny*nz);
    sol_1n.clear(); sol_1n.set_local(0, 0.1);

    Icarus::assembleFem assembler_1n(h, nx, ny, nz);
    assembler_1n.assemble(matrix_1n, rhs_1n, f_1, d_1, n_1);

    Icarus::BiCgStabSolver<Icarus::DistEllpackMatrix<double>> solver_1n(matrix_1n, rhs_1n);
    solver_1n.solve(sol_1n);

    Icarus::FullVector<double> fullsol_1n(sol_1n);
    Icarus::vtkWriter writer_1n("u=(0.5-x)...", "u=(0.5-x)...", nx, ny, nz, h, 1);
    writer_1n.addPointDataToTimestep(fullsol_1n, 0, "Potential");

    double l2norm_1n(0.0);
    for (int x(0); x < nx; ++x)
        for (int y(0); y < ny; ++y)
            for (int z(0); z < nz; ++z)
                l2norm_1n += (u_1.eval(x*h, y*h, z*h)-fullsol_1n[z*ny*nx+y*nx+x])*(u_1.eval(x*h, y*h, z*h)-fullsol_1n[z*ny*nx+y*nx+x]);
    l2norm_1n = sqrt(l2norm_1n);
    LOG_INFO("l2norm_1n = ", l2norm_1n);
*/
/*
    // ***** ***** ***** ***** TEST 2 DIRICHLET ***** ***** ***** ***** //
    // u = x*(1.0-x)*y*(1.0-y)*z*(1.0-z)
    // f = 2.0*(x*(1.0-x)*y*(1.0-y)+x*(1.0-x)*z*(1.0-z)+y*(1.0-y)*z*(1.0-z))
    // d = 0
    Icarus::mathfunction u_2(1);
    Icarus::mathfunction f_2(2);
    Icarus::mathfunction g_2(3);

    Icarus::DistEllpackMatrix<double> matrix_2d(nx*ny*nz);
    Icarus::SlicedVector<double> rhs_2d(nx*ny*nz);
    Icarus::SlicedVector<double> sol_2d(nx*ny*nz);
    sol_2d.clear(); sol_2d.set_local(0, 0.1);

    Icarus::assembleFem assembler_2d(h, nx, ny, nz);
    assembler_2d.assemble(matrix_2d, rhs_2d, f_2, g_2);

    Icarus::BiCgStabSolver<Icarus::DistEllpackMatrix<double>> solver_2d(matrix_2d, rhs_2d);
    solver_2d.solve(sol_2d);

    Icarus::FullVector<double> fullsol_2d(sol_2d);
    Icarus::vtkWriter writer_2d("u=(0.5-x)...", "u=(0.5-x)...", nx, ny, nz, h, 1);
    writer_2d.addPointDataToTimestep(fullsol_2d, 0, "Potential");

    double l2norm_2d(0.0);
    for (int x(0); x < nx; ++x)
        for (int y(0); y < ny; ++y)
            for (int z(0); z < nz; ++z)
                l2norm_2d += (u_2.eval(x*h, y*h, z*h)-fullsol_2d[z*ny*nx+y*nx+x])*(u_2.eval(x*h, y*h, z*h)-fullsol_2d[z*ny*nx+y*nx+x]);
    l2norm_2d = sqrt(l2norm_2d);
    LOG_INFO("l2norm_2d = ", l2norm_2d);

    // ***** ***** ***** ***** TEST 2 DIRICHLET ***** ***** ***** ***** //
    // u = x*(1.0-x)*y*(1.0-y)*z*(1.0-z)
    // f = 2.0*(x*(1.0-x)*y*(1.0-y)+x*(1.0-x)*z*(1.0-z)+y*(1.0-y)*z*(1.0-z))
    // d = 0
    // n = Fallunterscheidung

    LOG_INFO("l2norm_0d = ", l2norm_0d);
    LOG_INFO("l2norm_0n = ", l2norm_0n);
    LOG_INFO("l2norm_1d = ", l2norm_1d);
    LOG_INFO("l2norm_1n = ", l2norm_1n);
    LOG_INFO("l2norm_2d = ", l2norm_2d);
*/
    return 0;
}
