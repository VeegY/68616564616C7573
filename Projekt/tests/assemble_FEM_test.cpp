#include "../src/include/assemblefem.hpp"
#include "../src/include/mathfunction.hpp"

#include "../src/include/fullvector.hpp"
#include "../src/include/slicedvector.hpp"
#include "../src/include/bicgstabsolver.hpp"
#include "../src/include/distellpackmatrix.hpp"

#include "../src/include/vtkwriter.hpp"
#include "../src/include/discretizer.hpp"

#include <string>

void test(int testcase, bool neumann, std::vector<char>& disc, double h, int nx, int ny, int nz);

int main()
{
    const int nn(10);
    const double h(1.0/static_cast<double>(nn));
    const int nx(nn+1), ny(nn+1), nz(nn+1);
    std::vector<char> disc = Icarus::discretizer("../model/cube.obj", h, nx, ny, nz);

//    for (int i(1); i <= 3; ++i)
    {
        test(3, true, disc, h, nx, ny, nz);
    }

/*
    // ***** ***** ***** ***** TEST 0 DIRICHLET ***** ***** ***** ***** //
    // u = 0
    // f = 0
    // d = 0
    Icarus::mathfunction u_0(0, 0.0);
    Icarus::mathfunction f_0(0, 0.0);
    Icarus::mathfunction d_0(0, 0.0);

    Icarus::DistEllpackMatrix<double> matrix_0d(nx*ny*nz);
    Icarus::SlicedVector<double> rhs_0d(nx*ny*nz);
    Icarus::SlicedVector<double> sol_0d(nx*ny*nz);
    sol_0d.clear(); sol_0d.set_local(0, 0.1);

    Icarus::assembleFem assembler_0d(h, nx, ny, nz);
    assembler_0d.assemble(matrix_0d, rhs_0d, disc, f_0, d_0);

    Icarus::BiCgStabSolver<Icarus::DistEllpackMatrix<double>> solver_0d(matrix_0d, rhs_0d);
    solver_0d.solve(sol_0d);

    Icarus::FullVector<double> fullsol_0d(sol_0d);
//    Icarus::vtkWriter writer_0d("case0d", "case0d", nx+1, ny+1, nz+1, h, 1);
    Icarus::vtkWriter writer_0d("case0d", "case0d", nx, ny, nz, h, 1);
    writer_0d.addPointDataToTimestep(fullsol_0d, 0, "Potential");
//    writer_0d.addCellDataToTimestep(fullsol_0d, 0, "Potential");

    LOG_INFO("L2norm_0d = ", assembler_0d.calcL2Error(u_0, fullsol_0d));

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
    assembler_0n.assemble(matrix_0n, rhs_0n, disc, f_0, d_0, n_0);

    Icarus::BiCgStabSolver<Icarus::DistEllpackMatrix<double>> solver_0n(matrix_0n, rhs_0n);
    solver_0n.solve(sol_0n);

    Icarus::FullVector<double> fullsol_0n(sol_0n);
    Icarus::vtkWriter writer_0n("case0n", "case0n", nx, ny, nz, h, 1);
    writer_0n.addPointDataToTimestep(fullsol_0n, 0, "Potential");

    LOG_INFO("L2norm_0n = ", assembler_0n.calcL2Error(u_0, fullsol_0n));
*/
    return 0;
}

void test(int testcase, bool neumann, std::vector<char>& disc, double h, int nx, int ny, int nz)
{
    Icarus::mathfunction u((testcase-1)*4+1);
    Icarus::mathfunction f((testcase-1)*4+2);
    Icarus::mathfunction d((testcase-1)*4+3);
    Icarus::mathfunction n((testcase-1)*4+4);

    Icarus::DistEllpackMatrix<double> matrix(nx*ny*nz);
    Icarus::SlicedVector<double> rhs(nx*ny*nz);
    Icarus::SlicedVector<double> sol(nx*ny*nz);
    sol.clear(); sol.set_local(0, 0.1);

    Icarus::assembleFem assembler(h, nx, ny, nz);
    assembler.assemble(matrix, rhs, disc, f, d, n);

    Icarus::BiCgStabSolver<Icarus::DistEllpackMatrix<double>> solver(matrix, rhs);
    solver.solve(sol);

    Icarus::FullVector<double> fullsol(sol);
    std::string casename = "case" + std::to_string(testcase) + (neumann ? "n" : "d");
    Icarus::vtkWriter writer(casename, casename, nx, ny, nz, h, 1);
    writer.addPointDataToTimestep(fullsol, 0, "Potential");

    LOG_INFO("L2norm_", testcase, (neumann ? "n" : "d"), " = ", assembler.calcL2Error(u, fullsol));
}
