#include "../src/assemblyLGS.cpp"

#include "../src/include/slicedvector.hpp"
#include "../src/include/bicgstabsolver.hpp"
#include "../src/include/distellpackmatrix.hpp"

//#include "../src/include/vtkwriter.hpp"
//#include "../src/include/discretizer.hpp"

int main()
{
    Icarus::DistEllpackMatrix<double> matrix(Nx*Ny*Nz);

    Icarus::SlicedVector<double> rhs(Nx*Ny*Nz);
    rhs.clear();

    Icarus::SlicedVector<double> sol(Nx*Ny*Nz);
    sol.clear();

    Icarus::assemble_FEM(matrix, rhs);
/*
	Icarus::BiCgStabSolver<Icarus::DistEllpackMatrix<double>> solver(matrix, rhs);
    solver.solve(sol);
    */
    return 0;
}
