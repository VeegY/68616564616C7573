#include "../src/include/fullvector.hpp"
#include "../src/include/vtkwriter.hpp"
/*
*test fuer den vtk writer
*ausgabe im ordner out (muss existieren)
*KEIN unit-test
*
*/
int gradient_test()
{
	const int Nx = 5, Ny = 5, Nz = 5, points=Nx*Ny*Nz, cells=(Nx-1)*(Ny-1)*(Nz-1);
	const float h = 0.3;

	Icarus::FullVector<double> potential(points), pointx(points), pointy(points), pointz(points);
	Icarus::FullVector<float> cellx(cells), celly(cells), cellz(cells);
	double h = 0.2;
	size_t globIdx(0);
	for (size_t z(0); z < nz; z++)
    {
        for (size_t y(0); y < ny; y++)
        {
            for (size_t x(0); x < nx; x++ , globIdx++)
            {
                potential=y*z*x*h*h*h;
            }
        }
    }
    Icarus::getInnerPotentialGradients(potential, Nx, Ny, Nz, h, pointx, pointy, pointz);
    Icarus::getCellMidpointGradientsFEM(potential, Nx, Ny, Nz, cellx, celly, cellz);
    Icarus::vtkWriter writer("gradienttest", "test", Nx, Ny, Nz, h, 1);
    writer.addPointDataToAll(potential, "potential");
    writer.addPointVecToTimestep(pointx, pointy, pointz,0, "inner gradients");
    writer.addCellVecToTimestep(cellx, celly, cellz, 0, "cell midpoint gradients");
    return 0;
}

int main()
{
    int myrank;
    MPI_SCALL(MPI_Comm_rank(MPI_COMM_WORLD, &myrank));
    if (myrank == 0)
    {
        return gradient_test();
    }
    return 0;
}

