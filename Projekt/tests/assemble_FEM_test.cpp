#include "../src/include/assemblefem.hpp"
#include "../src/include/mathfunction.hpp"

#include "../src/include/fullvector.hpp"
#include "../src/include/slicedvector.hpp"
#include "../src/include/bicgstabsolver.hpp"
#include "../src/include/distellpackmatrix.hpp"

#include "../src/include/vtkwriter.hpp"
#include "../src/include/discretizer.hpp"

int main()
{
    const int nn(10);
    const double h(1.0/static_cast<double>(nn));
    const int nx(nn+1), ny(nn+1), nz(nn+1);
    std::vector<char> disc = Icarus::discretizer("leer.obj", h, nx, ny, nz);

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
    assembler_0d.assemble(matrix_0d, rhs_0d, disc, f_0, d_0);

    Icarus::BiCgStabSolver<Icarus::DistEllpackMatrix<double>> solver_0d(matrix_0d, rhs_0d);
    solver_0d.solve(sol_0d);

    Icarus::FullVector<double> fullsol_0d(sol_0d);
    Icarus::vtkWriter writer_0d("case0d", "case0d", nx, ny, nz, h, 1);
    writer_0d.addPointDataToTimestep(fullsol_0d, 0, "Potential");

    LOG_INFO("L2norm_0d = ", assembler_0d.calcL2Error(u_0, fullsol_0d));
*/
    // ***** ***** ***** ***** TEST 0 NEUMANN ***** ***** ***** ***** //
    // u = 0
    // f = 0
    // d = 0
    // n = 0
    Icarus::mathfunction n_0(0, 0.0);
/*
    Icarus::DistEllpackMatrix<double> matrix_0n(nx*ny*nz);
    Icarus::SlicedVector<double> rhs_0n(nx*ny*nz);
    Icarus::SlicedVector<double> sol_0n(nx*ny*nz);
    sol_0n.clear(); sol_0n.set_local(0, 0.1);

    Icarus::assembleFem assembler_0n(h, nx, ny, nz);
    assembler_0n.assemble(matrix_0n, rhs_0n, disc, f_0, d_0, n_0);

    Icarus::BiCgStabSolver<Icarus::DistEllpackMatrix<double>> solver_0n(matrix_0n, rhs_0n);
    solver_0n.solve(sol_0n);

    Icarus::FullVector<double> fullsol_0n(sol_0n);
    Icarus::vtkWriter writer_0n("case0n", "case0n", nx, ny, nz, h, 1);
    writer_0n.addPointDataToTimestep(fullsol_0n, 0, "Potential");

    LOG_INFO("L2norm_0n = ", assembler_0n.calcL2Error(u_0, fullsol_0n));
*/
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
    assembler_1d.assemble(matrix_1d, rhs_1d, disc, f_1, d_1);

    Icarus::BiCgStabSolver<Icarus::DistEllpackMatrix<double>> solver_1d(matrix_1d, rhs_1d);
    solver_1d.solve(sol_1d);

    Icarus::FullVector<double> fullsol_1d(sol_1d);
    Icarus::vtkWriter writer_1d("case1d", "case1d", nx, ny, nz, h, 1);
    writer_1d.addPointDataToTimestep(fullsol_1d, 0, "Potential");

    LOG_INFO("L2norm_1d = ", assembler_1d.calcL2Error(u_1, fullsol_1d));
*/
    // ***** ***** ***** ***** TEST 1 NEUMANN ***** ***** ***** ***** //
    // u = (0.5-x)(0.5-y)(0.5-z)
    // f = 0
    // d = Fallunterscheidung...
    // n = Fallunterscheidung...
    Icarus::mathfunction n_1(4);

    Icarus::DistEllpackMatrix<double> matrix_1n(nx*ny*nz);
    Icarus::SlicedVector<double> rhs_1n(nx*ny*nz);
    Icarus::SlicedVector<double> sol_1n(nx*ny*nz);
    sol_1n.clear(); sol_1n.set_local(0, 0.1);

    Icarus::assembleFem assembler_1n(h, nx, ny, nz);
    assembler_1n.assemble(matrix_1n, rhs_1n, disc, f_1, d_1, n_1);

    Icarus::BiCgStabSolver<Icarus::DistEllpackMatrix<double>> solver_1n(matrix_1n, rhs_1n);
    solver_1n.solve(sol_1n);

    Icarus::FullVector<double> fullsol_1n(sol_1n);
    Icarus::vtkWriter writer_1n("case1n", "case1n", nx, ny, nz, h, 1);
    writer_1n.addPointDataToTimestep(fullsol_1n, 0, "Potential");

    LOG_INFO("L2norm_1n = ", assembler_1n.calcL2Error(u_1, fullsol_1n));
//    matrix_1n.print_local_data(std::cout);

    // ***** ***** ***** ***** TEST 2 DIRICHLET ***** ***** ***** ***** //
    // u = x*(1.0-x)*y*(1.0-y)*z*(1.0-z)
    // f = 2.0*(x*(1.0-x)*y*(1.0-y)+x*(1.0-x)*z*(1.0-z)+y*(1.0-y)*z*(1.0-z))
    // d = 0
    Icarus::mathfunction u_2(5);
    Icarus::mathfunction f_2(6);
    Icarus::mathfunction d_2(7);
/*
    Icarus::DistEllpackMatrix<double> matrix_2d(nx*ny*nz);
    Icarus::SlicedVector<double> rhs_2d(nx*ny*nz);
    Icarus::SlicedVector<double> sol_2d(nx*ny*nz);
    sol_2d.clear(); sol_2d.set_local(0, 0.1);

    Icarus::assembleFem assembler_2d(h, nx, ny, nz);
    assembler_2d.assemble(matrix_2d, rhs_2d, disc, f_2, d_2);

    Icarus::BiCgStabSolver<Icarus::DistEllpackMatrix<double>> solver_2d(matrix_2d, rhs_2d);
    solver_2d.solve(sol_2d);

    Icarus::FullVector<double> fullsol_2d(sol_2d);
    Icarus::vtkWriter writer_2d("case2d", "case2d", nx, ny, nz, h, 1);
    writer_2d.addPointDataToTimestep(fullsol_2d, 0, "Potential");

    LOG_INFO("L2norm_2d = ", assembler_2d.calcL2Error(u_2, fullsol_2d));
*/
    // ***** ***** ***** ***** TEST 2 NEUMANN ***** ***** ***** ***** //
    // u = x*(1.0-x)*y*(1.0-y)*z*(1.0-z)
    // f = 2.0*(x*(1.0-x)*y*(1.0-y)+x*(1.0-x)*z*(1.0-z)+y*(1.0-y)*z*(1.0-z))
    // d = 0
    // n = Fallunterscheidung
    Icarus::mathfunction n_2(8);
/*
    Icarus::DistEllpackMatrix<double> matrix_2n(nx*ny*nz);
    Icarus::SlicedVector<double> rhs_2n(nx*ny*nz);
    Icarus::SlicedVector<double> sol_2n(nx*ny*nz);
    sol_2n.clear(); sol_2n.set_local(0, 0.1);

    Icarus::assembleFem assembler_2n(h, nx, ny, nz);
    assembler_2n.assemble(matrix_2n, rhs_2n, disc, f_2, d_2, n_2);

    Icarus::BiCgStabSolver<Icarus::DistEllpackMatrix<double>> solver_2n(matrix_2n, rhs_2n);
    solver_2n.solve(sol_2n);

    Icarus::FullVector<double> fullsol_2n(sol_2n);
    Icarus::vtkWriter writer_2n("case2n", "case2n", nx, ny, nz, h, 1);
    writer_2n.addPointDataToTimestep(fullsol_2n, 0, "Potential");

    LOG_INFO("L2norm_2n = ", assembler_2n.calcL2Error(u_2, fullsol_2n));
//    matrix_2n.print_local_data(std::cout);
*/
    return 0;
}
