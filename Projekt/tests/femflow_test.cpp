#include "../src/include/assemblefem.hpp"
#include "../src/include/mathfunction.hpp"

#include "../src/include/fullvector.hpp"
#include "../src/include/slicedvector.hpp"
#include "../src/include/bicgstabsolver.hpp"
#include "../src/include/distellpackmatrix.hpp"

#include "../src/include/fullvectorgpu.hpp"
#include "../src/include/slicedvectorgpu.hpp"
#include "../src/include/distellpackmatrixgpu.hpp"

#include "../src/include/vtkwriter.hpp"
#include "../src/include/discretizer.hpp"
#include "../src/include/potentialgradients.hpp"

int main()
{
    const int nn(50);
    const double h(1.0/static_cast<double>(nn));
    const int nx(nn+1), ny(nn+1), nz(nn+1);
    std::vector<char> disc = Icarus::discretizer("../model/cube.obj", h, nx, ny, nz);

    Icarus::mathfunction u(13);
    Icarus::mathfunction f(14);
    Icarus::mathfunction d(15);
    Icarus::mathfunction n(16);

    // CPU
    Icarus::DistEllpackMatrix<double> matrix(nx*ny*nz);
    Icarus::SlicedVector<double> rhs(nx*ny*nz);
    Icarus::SlicedVector<double> sol(nx*ny*nz);
    sol.clear(); sol.set_local(0, 0.1);
    Icarus::BiCgStabSolver<Icarus::DistEllpackMatrix<double>> solver(matrix, rhs);
    Icarus::assembleFem assembler(h, nx, ny, nz);

    assembler.assemble(matrix, rhs, disc, f, d, n);
    solver.solve(sol);

    Icarus::FullVector<double> fullsol(sol);
    Icarus::FullVector<double> gradx(nx*ny*nz), grady(nx*ny*nz), gradz(nx*ny*nz);
    Icarus::getInnerPotentialGradients(fullsol, nx, ny, nz, h, disc, gradx, grady, gradz);
    Icarus::vtkWriter writer("casefemflow", "casefemflow", nx, ny, nz, h, 1);
    writer.addPointDataToTimestep(fullsol, 0, "Potential");
    writer.addPointVecToTimestep(gradx, grady, gradz, 0, "Gradient");

/*
    // GPU
    Icarus::DistEllpackMatrixGpu<double> matrix(nx*ny*nz);
    Icarus::SlicedVectorGpu<double> rhs(nx*ny*nz);
    Icarus::SlicedVectorGpu<double> sol(nx*ny*nz);
    sol.clear(); sol.set_local(0, 0.1);
    Icarus::BiCgStabSolver<Icarus::DistEllpackMatrixGpu<double>> solver(matrix, rhs);
    Icarus::assembleFem assembler(h, nx, ny, nz);

    assembler.assemble(matrix, rhs, disc, f, d, n);
    solver.solve(sol);

    Icarus::FullVectorGpu<double> fullsol(sol);
    Icarus::FullVectorGpu<double> gradx(nx*ny*nz), grady(nx*ny*nz), gradz(nx*ny*nz);

    Icarus::getInnerPotentialGradients(fullsol, nx, ny, nz, h, disc, gradx, grady, gradz);
    Icarus::vtkWriter writer("casefemflow", "casefemflow", nx, ny, nz, h, 1);
    writer.addPointDataToTimestep(fullsol.getDataPointer(), fullsol.get_dim(), 0, "Potential");
    writer.addPointVecToTimestep(gradx.getDataPointer(), grady.getDataPointer(), gradz.getDataPointer(), gradx.get_dim(), 0, "Gradient");
*/

    return 0;
}
