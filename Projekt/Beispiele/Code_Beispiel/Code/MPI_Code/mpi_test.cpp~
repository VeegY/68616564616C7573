// #################################################################################################
//					 FD - Code Projekt
// #################################################################################################
// 				     CPP-File: mpi_test
// ------------------------------------Doxygen-Dokumentation----------------------------------------
///  \file mpi_test.cpp
///  \brief
///  Fuehrt eine 'Hallo Welt'-MPI Implementierung aus
// -------------------------------------------------------------------------------------------------
// Kompilieren mit dem folgenden Compilerbefehl: !! WIP !!
//   mpicc mpi_test.cpp
// -------------------------------------------------------------------------------------------------
// Verwendete Header:
#	include <iostream>
#	include <mpi>

//#	include ".h"
// -------------------------------------------------------------------------------------------------
// Copyright : 2015 malte-christian.loddoch@tu-dortmund.de
// #################################################################################################

static int numprocs;

// =================================================================================================
//					 main-Funktion
// -------------------------------------------------------------------------------------------------
///  \brief Fuehrt den 'Hallo Welt' MPI-Test durch.
// =================================================================================================
int main(int argc, char **argv)
{ 
	int my_rank; 
	MPI_Status status; 	// MPI initializations 
	MPI_Init (&argc, &argv); 
	MPI_Comm_size (MPI_COMM_WORLD, &numprocs); 
	MPI_Comm_rank (MPI_COMM_WORLD, &my_rank); 
	double time_start = MPI_Wtime();
	std::cout << "Hello World, my rank is " << my_rank <<" "<< MPI_Wtime() - time_start << std::endl; 
	MPI_Finalize (); 	// End MPI 
	return 0; 
}
