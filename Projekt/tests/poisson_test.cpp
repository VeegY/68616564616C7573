/*
#include "../src/include/slicedvector.hpp"
#include "../src/include/bicgstabsolver.hpp"
#include "../src/include/distellpackmatrix.hpp"

#include "../src/include/vtkwriter.hpp"
#include "../src/include/discretizer.hpp"
#include "../src/include/assemble.hpp"

double rhs(int ix, int iy, int iz)
{
	if (iz <= 20 && ix <= 20 && iy <= 20) return 20.0;
	else if (iz >= 30 && ix >= 30 && iy >= 30) return -20.0;
	return 0.0;
}

int poisson_demo()
{
	const int Nx = 50, Ny = 50, Nz = 50;
	const float h = 0.1;
	// diskretisieren
	std::vector<char> disc = Icarus::discretizer("leer.obj", h, Nx, Ny, Nz);
	
	// assemblieren
	auto lgs = Icarus::assemble<double>(disc, h, Nx, Ny, Nz, rhs);
	
	// lösen
	size_t n = lgs.first.get_dim_global();
	Icarus::SlicedVector<double> sol(n);
	sol.clear();
	Icarus::BiCgStabSolver<Icarus::DistEllpackMatrix<double>> solver(lgs.first, lgs.second);
	solver.solve(sol);

	// speichern (nur master)
	Icarus::FullVector<double> fullsol(sol);	
	int myrank;
	MPI_SCALL(MPI_Comm_rank(MPI_COMM_WORLD, &myrank));
	if (myrank == 0)
	{
		Icarus::vtkWriter writer("out", "Testdatensatz", Nx, Ny, Nz, 1);
		writer.addPointDataToTimestep(fullsol, 0, "Potential");
	}
	MPI_SCALL(MPI_Barrier(MPI_COMM_WORLD));
	return 0;
}

int main()
{
    return poisson_demo();
}
*/
int main()
{
    return 1;
}
