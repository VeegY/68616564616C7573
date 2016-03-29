
#include "../src/include/slicedvector.hpp"
#include "../src/include/bicgstabsolver.hpp"
#include "../src/include/distellpackmatrix.hpp"

#include "../src/include/vtkwriter.hpp"
#include "../src/include/discretizer.hpp"
#include "../src/include/assemble.hpp"

double bdry(int vtx_global)
{
	if (vtx_global < 100) return 10000;
	if (vtx_global >= 900) return -10000;
	return 0.0;
}

int neumann_demo()
{
	const int nx = 20, ny = 20, nz = 20;
	const float h = 0.05;
	// diskretisieren
	//std::vector<char> disc = Icarus::discretizer("leer.obj", h, Nx, Ny, Nz);
	
	// assemblieren
	auto lgs = Icarus::assemble_neumann_unrolled<double>(nx, ny, nz, h, bdry);
	//ausgabe der Matrix:
	//std::cout << "Die Matrix A:" << std::endl;
	//Icarus::print_sliced_object(lgs.first);
	
	// loesen
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
		Icarus::vtkWriter writer("out/neumann", "Testdatensatz", nx, ny, nz, 1);
		writer.addPointDataToTimestep(fullsol, 0, "Geschwindigkeitspotential");
	}
	MPI_SCALL(MPI_Barrier(MPI_COMM_WORLD));
	return 0;
}

int main()
{
	return neumann_demo();
}

