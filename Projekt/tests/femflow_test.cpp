#include "../src/include/assemblefem.hpp"
#include "../src/include/mathfunction.hpp"

#include "../src/include/fullvector.hpp"
#include "../src/include/slicedvector.hpp"
#include "../src/include/bicgstabsolver.hpp"
#include "../src/include/distellpackmatrix.hpp"

#include "../src/include/vtkwriter.hpp"
#include "../src/include/discretizer.hpp"
#include "../src/include/potentialgradients.hpp"

int main()
{
    const int nn(100);
    const double h(1.0/static_cast<double>(nn));
    const int nx(nn+1), ny(nn+1), nz(nn+1);
    std::vector<char> disc = Icarus::discretizer("../model/board_big.obj", h, nx, ny, nz);

    Icarus::mathfunction u(13);
    Icarus::mathfunction f(14);
    Icarus::mathfunction d(15);
    Icarus::mathfunction n(16);

    Icarus::DistEllpackMatrix<double> matrix(nx*ny*nz);
    Icarus::SlicedVector<double> rhs(nx*ny*nz);
    Icarus::SlicedVector<double> sol(nx*ny*nz);
    sol.clear(); sol.set_local(0, 0.1);

    Icarus::assembleFem assembler(h, nx, ny, nz);
    assembler.assemble(matrix, rhs, disc, f, d, n);

    Icarus::BiCgStabSolver<Icarus::DistEllpackMatrix<double>> solver(matrix, rhs);
    solver.solve(sol);

    Icarus::FullVector<double> fullsol(sol);
    Icarus::FullVector<double> gradx(nx*ny*nz), grady(nx*ny*nz), gradz(nx*ny*nz);
    Icarus::getInnerPotentialGradients(fullsol, nx, ny, nz, h, disc, gradx, grady, gradz);
    Icarus::vtkWriter writer("casefemflow", "casefemflow", nx, ny, nz, h, 1);
    writer.addPointDataToTimestep(fullsol, 0, "Potential");
    writer.addPointVecToTimestep(gradx, grady, gradz, 0, "Gradient");

//    LOG_INFO("L2norm = ", assembler.calcL2Error(u, fullsol));

    return 0;
}
