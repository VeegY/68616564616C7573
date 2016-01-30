#include "include\fullvector.hpp"
#include "include\vtkwriter.hpp"
/*
*test für den vtk writer
*ausgabe im ordner out (muss existieren)
*KEIN unit-test
*
*/
int vtk_writer_test()
{
	const int Nx = 3, Ny = 3, Nz = 3, points=27, cells=8;
	const float h = 0.3;

	Icarus::FullVector<double> pointx(27), pointy(27),pointz(27);
	Icarus::FullVector<float> cellx(8), celly(8), cellyz(8);
	float arr_pointx[27], arr_pointy[27], arr_pointz[27];
	double arr_cellx[8], arr_celly[8], arr_cellz[8];

	//set values
	for (int i(0); i<points; i++){
        pointx[i]=i;
        pointy[i]=i*i;
        pointz[i]=27-i;
        arr_pointx[i]=i;
        arr_pointy[i]=i*i;
        arr_pointz[i]=27-i;
	}

	for (int i(0); i<8; i++){
        cellx[i]=i;
        celly[i]=i*i;
        cellz[i]=8-i;
        arr_cellx[i]=i;
        arr_celly[i]=i*i;
        arr_cellz[i]=8-i;
	}
    Icarus::vtkWriter writer("out/writer_test", "Testdatensatz", Nx, Ny, Nz, 1);
    Icarus::vtkWriter writer2("out/writer_test2", "Testdatensatz", Nx, Ny, Nz, h, 1);

    writer.addPointDataToAll(pointx, "point");
    writer2.addPointDataToAll(arrpointx, points, "point");

    writer.addCellDataToAll(cellx, "cell");
    writer2.addCellDataToAll(arrcellx, cells, "cell");

    writer.addPointVecToAll(pointx, pointy, pointz, "pointvec");
    writer2.addPointVecToAll(arrpointx, arrpointy, arrpointz, points, "pointvec");

    writer.addCellVecToAll(cellx, celly, cellz, "cellvec");
    writer2.addCellVecToAll(arrcellx, arrcelly, arrcelz, cells, "cellvec");

	return 0;
}

