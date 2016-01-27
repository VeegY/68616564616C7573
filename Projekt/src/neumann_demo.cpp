#include "include\slicedvector.hpp"
#include "include\bicgstabsolver.hpp"
#include "include\distellpackmatrix.hpp"

#include "include\vtkwriter.hpp"
#include "include\discretizer.hpp"
#include "include\assemble.hpp"

double bdry(int vtx_global)
{
	if (vtx_global < 100) return 100;
	return 0.0;
}

int neumann_demo()
{
	const int nx = 50, ny = 50, nz = 50;
	const float h = 0.1;
	// diskretisieren
	//std::vector<char> disc = Icarus::discretizer("leer.obj", h, Nx, Ny, Nz);
	
	// assemblieren
	auto lgs = Icarus::assemble_neumann<double>(nx, ny, nz, h, bdry);
	
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
		Icarus::vtkWriter writer("out", "Testdatensatz", nx, ny, nz, 1);
		writer.addPointDataToTimestep(fullsol, 0, "Geschwindigkeitspotential");
	}
	MPI_SCALL(MPI_Barrier(MPI_COMM_WORLD));
	return 0;
}

int main()
{
	return neumann_demo();
}