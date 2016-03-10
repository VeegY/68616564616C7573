#include "../src/include/slicedvector.hpp"
#include "../src/include/bicgstabsolver.hpp"
#include "../src/include/distellpackmatrix.hpp"

#include "../src/include/logger.hpp"

#include <cmath>

int main()
{
    double TOL(1.0e-8);
    Icarus::DistEllpackMatrix<double> matrix(3);
    matrix.prepare_sequential_fill(3);
    matrix.sequential_fill(0, 1.0);
    matrix.end_of_row();
    matrix.sequential_fill(0, -1.0);
    matrix.sequential_fill(1, 10.0);
    matrix.sequential_fill(2, -1.0);
    matrix.end_of_row();
    matrix.sequential_fill(2, 1.0);
    matrix.end_of_row();

    //TODO Loeser schmeisst immer einen Error. Aber hier soll ja eigentlich auch nichts konvergieren. Kann man das umgehen,
    // sodass das hier failt, wenn es konvergiert und nicht failt, wenn es nicht konvergiert?
    // Anstonsten bleibt das noch zum per Hand Testen auskommentiert drin.
    // Rechte Seite, die in Matlab nicht konvergiert
//    Icarus::SlicedVector<double> rhs1(3);
//    rhs1.set_global(0, 0.0);
//    rhs1.set_global(1, 1.0);
//    rhs1.set_global(2, 0.0);
//
//    Icarus::SlicedVector<double> sol1(3);
//    Icarus::BiCgStabSolver<Icarus::DistEllpackMatrix<double>> solver1(matrix, rhs1);
//    solver1.solve(sol1);
//
//    Icarus::SlicedVector<double> correctsol1(3);
//    correctsol1.set_global(0, 0.0);
//    correctsol1.set_global(1, 0.1);
//    correctsol1.set_global(2, 0.0);
//
//    //TODO direkt was mit Norm?
//    double l2norm1(sqrt((sol1.get_global(0) - correctsol1.get_global(0)) * (sol1.get_global(0) - correctsol1.get_global(0)) +
//                        (sol1.get_global(1) - correctsol1.get_global(1)) * (sol1.get_global(1) - correctsol1.get_global(1)) +
//                        (sol1.get_global(2) - correctsol1.get_global(2)) * (sol1.get_global(2) - correctsol1.get_global(2))));
//    if (l2norm1 < TOL)
//        LOG_INFO("Konvergenz zur eindeutigen Loesung, obwohl Matlabs BiCGStab nicht konvergiert.");
//    else
//        LOG_INFO("Keine Konvergenz zur eindeutigen Loesung. Matlab BiCGStab konvergiert in dem Fall aber auch nicht.");

    // Rechte Seite, die in Matlab konvergiert
    Icarus::SlicedVector<double> rhs2(3);
    rhs2.set_global(0, 1.0);
    rhs2.set_global(1, 0.0);
    rhs2.set_global(2, 1.0);
    Icarus::print_sliced_object(rhs2);
    Icarus::SlicedVector<double> sol2(3);
    Icarus::BiCgStabSolver<Icarus::DistEllpackMatrix<double>> solver2(matrix, rhs2);
    sol2.clear();
    sol2.set_local(0, 1e-6);
    solver2.solve(sol2);

    Icarus::SlicedVector<double> correctsol2(3);
    correctsol2.set_global(0, 1.0);
    correctsol2.set_global(1, 0.2);
    correctsol2.set_global(2, 1.0);

    //TODO direkt was mit Norm?
    double l2norm2(sqrt((sol2.get_global(0) - correctsol2.get_global(0)) * (sol2.get_global(0) - correctsol2.get_global(0)) +
                        (sol2.get_global(1) - correctsol2.get_global(1)) * (sol2.get_global(1) - correctsol2.get_global(1)) +
                        (sol2.get_global(2) - correctsol2.get_global(2)) * (sol2.get_global(2) - correctsol2.get_global(2))));
    if (l2norm2 < TOL)
        LOG_INFO("Konvergenz zur eindeutigen Loesung. Matlab BiCGStab konvergiert in dem Fall auch.");
    else
        LOG_ERROR("Keine Konvergenz zur eindeutigen Loesung, obwohl Matlabs BiCGStab konvergiert.");

    return 0;
}
