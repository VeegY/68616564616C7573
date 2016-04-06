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
    LOG_INFO("start CPU test");
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
    Icarus::vtkWriter writer("casefemflow_cpu", "casefemflow_cpu", nx, ny, nz, h, 1);
    writer.addPointDataToTimestep(fullsol, 0, "Potential");
    writer.addPointVecToTimestep(gradx, grady, gradz, 0, "Gradient");
    LOG_INFO("end CPU test");

    // GPU
    LOG_INFO("start GPU test");
    Icarus::DistEllpackMatrixGpu<double> matrix_gpu(nx*ny*nz);
    Icarus::SlicedVectorGpu<double> rhs_gpu(nx*ny*nz);
    Icarus::SlicedVectorGpu<double> sol_gpu(nx*ny*nz);
    sol_gpu.clear(); sol_gpu.set_local(0, 0.1);
    Icarus::BiCgStabSolver<Icarus::DistEllpackMatrixGpu<double>> solver_gpu(matrix_gpu, rhs_gpu);
    Icarus::assembleFem assembler_gpu(h, nx, ny, nz);

    assembler_gpu.assemble(matrix_gpu, rhs_gpu, disc, f, d, n);
    solver_gpu.solve(sol_gpu);

    Icarus::FullVectorGpu<double> fullsol_gpu(sol_gpu);
    Icarus::FullVectorGpu<double> gradx_gpu(nx*ny*nz), grady_gpu(nx*ny*nz), gradz_gpu(nx*ny*nz);

    Icarus::getInnerPotentialGradients(fullsol_gpu, nx, ny, nz, h, disc, gradx_gpu, grady_gpu, gradz_gpu);
    Icarus::vtkWriter writer_gpu("casefemflow_gpu", "casefemflow_gpu", nx, ny, nz, h, 1);
    writer.addPointDataToTimestep(fullsol_gpu.getDataPointer(), fullsol_gpu.get_dim(), 0, "Potential");
    writer.addPointVecToTimestep(gradx_gpu.getDataPointer(), grady_gpu.getDataPointer(), gradz_gpu.getDataPointer(), gradx_gpu.get_dim(), 0, "Gradient");
    LOG_INFO("end GPU test");

    return 0;
}
