#include "../src/include/fullvector.hpp"
#include "../src/include/vtkwriter.hpp"
/*
*test fuer den vtk writer
*ausgabe im ordner out (muss existieren)
*KEIN unit-test
*
*/
int vtk_writer_test()
{
	const int Nx = 3, Ny = 3, Nz = 3, points=27, cells=8;
	const float h = 0.3;

	Icarus::FullVector<double> pointx(27), pointy(27),pointz(27);
	Icarus::FullVector<float> cellx(8), celly(8), cellz(8);
	float arrpointx[27], arrpointy[27], arrpointz[27];
	double arrcellx[8], arrcelly[8], arrcellz[8];

	//set values
	for (int i(0); i<points; i++){
        pointx[i]=i;
        pointy[i]=i*i;
        pointz[i]=27-i;
        arrpointx[i]=i;
        arrpointy[i]=i*i;
        arrpointz[i]=27-i;
	}

	for (int i(0); i<8; i++){
        cellx[i]=i;
        celly[i]=i*i;
        cellz[i]=8-i;
        arrcellx[i]=i;
        arrcelly[i]=i*i;
        arrcellz[i]=8-i;
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
    writer2.addCellVecToAll(arrcellx, arrcelly, arrcellz, cells, "cellvec");

	return 0;
}

