#include "include\slicedvector.hpp"
#include "include\bicgstabsolver.hpp"
#include "include\distellpackmatrix.hpp"

#include "include\vtkwriter.hpp"
#include "include\discretizer.hpp"
#include "include\assemble.hpp"

double rhs(int ix, int iy, int iz)
{
	return 5.0 * (ix < 10 && iy < 10 && iz < 10);
}

int main()
{
	// diskretisieren
	std::vector<char> disc = Icarus::discretizer("leer.obj", 0.1, 50, 50, 50);
	
	// assemblieren
	auto lgs = Icarus::assemble<double>(disc, 0.1, 50, 50, 50, rhs);
	MPI_SCALL(MPI_Barrier(MPI_COMM_WORLD));
	
	// lösen
	size_t n = lgs.first.get_dim_global();
	Icarus::SlicedVector<double> sol(n);
	sol.clear();
	Icarus::BiCgStabSolver<Icarus::DistEllpackMatrix<double>> solver(lgs.first, lgs.second);
	solver.solve(sol);
	
	// speichern
	Icarus::FullVector<double> fullsol(sol);
	Icarus::vtkWriter writer("out", "Testdatensatz", 50, 50, 50, 1);
	writer.addPointDataToTimestep(fullsol, 0, "Potential");
	return 0;
}