#include "include/assemblefem.hpp"
#include "include/mathfunction.hpp"

#include "include/fullvector.hpp"
#include "include/slicedvector.hpp"
#include "include/bicgstabsolver.hpp"
#include "include/distellpackmatrix.hpp"

#include "include/fullvectorgpu.hpp"
#include "include/slicedvectorgpu.hpp"
#include "include/distellpackmatrixgpu.hpp"

#include "include/vtkwriter.hpp"
#include "include/discretizer.hpp"
#include "include/potentialgradients.hpp"

#include <string>

// Argumentuebergabe:
// Argumente in Klammern muessen nicht gesetzt werden.
// Wenn sie benutzt werden, muessen auch alle vorherigen gesetzt werden.
// ./useful_demo objfile stepsize (sizex) (sizey) (sizez) (gpu)
// objfile:  Dateiname ohne .obj Endung, welche im Ordner "model" liegen muss
//           string
// stepsize: Gibt die Schrittweite an
//           double
// sizex:    Ausdehnung in x-Richtung vom Ursprung aus (nicht Anzahl der Elemente!)
//           int/float/double, default: 1
// sizey:    Ausdehnung in y-Richtung vom Ursprung aus (nicht Anzahl der Elemente!)
//           int/float/double, default: 1
// sizez:    Ausdehnung in z-Richtung vom Ursprung aus (nicht Anzahl der Elemente!)
//           int/float/double, default: 1
// gpu:      Gibt an, ob auf der GPU oder auf der CPU gerechnet werden soll
//           "cpu" oder "gpu", default: cpu, TODO/ Bemerkung: gpu noch nicht eingebaut
// Beispiel:
// ./demo board 0.005 10 10 5 gpu
// TODO: Funktionennummern uebergeben koennen
// TODO: gpu Berechnung

int main(int argc, char* argv[])
{
    // setup durch Uebergabeparameter
    LOG_INFO("running program: ", argv[0]);
    if (argc < 3)
        LOG_ERROR("Zu wenig Argumente uebergeben.");
    std::string execpath(argv[0]);
    execpath.erase(0, 2);
    if (execpath.find_last_of('/') != std::string::npos)
        execpath.erase(execpath.find_last_of('/') + 1, execpath.size());
    else
        execpath = "";
    std::string model(argv[1]);
    std::string modelpath(execpath + "../../model/" + model + ".obj");
    LOG_INFO("use obj file: ", modelpath);
    //model += "_" + std::to_string(nn);
    const double h(std::stod(argv[2]));
    LOG_INFO("using step width: ", h);
    double sizex(1.0), sizey(1.0), sizez(1.0);
    if (argc > 3)
        sizex = std::stod(argv[3]);
    if (argc > 4)
        sizey = std::stod(argv[4]);
    if (argc > 5)
        sizez = std::stod(argv[5]);
    LOG_INFO("size: ", sizex, " x ", sizey, " x ", sizez);
    const int nx(sizex / h + 1), ny(sizey / h + 1), nz(sizez / h + 1);
    LOG_INFO("elements: ", nx-1, " x ", ny-1, " x ", nz-1); // Anzahl Elemente = Anzahl Stuetzstellen - 1
    LOG_INFO("calculate on cpu"); // TODO gpu ermoeglichen

    // mpi setup
    int myrank;
    MPI_SCALL(MPI_Comm_rank(MPI_COMM_WORLD, &myrank));

    std::vector<char> disc = Icarus::discretizer(modelpath, h, nx, ny, nz);

    // TODO: Funktionennummern uebergeben koennen
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

//matrix.~DistEllpackMatrix(); // um groessere Probleme berechnen zu koennen
    Icarus::FullVector<double> fullsol(sol);
    if (myrank == 0)
    {
        Icarus::FullVector<double> gradx(nx*ny*nz), grady(nx*ny*nz), gradz(nx*ny*nz);
        Icarus::getInnerPotentialGradients(fullsol, nx, ny, nz, h, disc, gradx, grady, gradz);
        Icarus::vtkWriter writer(model + "_cpu", model + "inner_cpu", nx, ny, nz, h, 1);
        writer.addPointDataToTimestep(fullsol, 0, "Potential");
        writer.addPointVecToTimestep(gradx, grady, gradz, 0, "Gradient");
        LOG_INFO("end CPU test");
    }
    MPI_SCALL(MPI_Barrier(MPI_COMM_WORLD));

    return 0;
}
