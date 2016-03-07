#ifndef __SOLVERBENCH_HPP_
#define __SOLVERBENCH_HPP_

#include <iostream>
#include <vector>
#include <mpi.h>
#include <chrono>
#include "bicgstabsolver.hpp"
#include "distellpackmatrix.hpp"

namespace Icarus
{
	// ACHTUNG: Nur die Prozesse, die in comm enthalten sind, dürfen diese Funktion aufrufen.
	template <class DistMatrix>
	static DistMatrix construct_model_matrix(unsigned m, MPI_Comm comm = MPI_COMM_WORLD) const
	{
		if (comm == MPI_COMM_NULL)
		{
			std::cerr "construct_model_matrix von ungueltigem comm aus aufgerufen!" << std::endl;
			throw;
		}

		DistMatrix mat(m*m, comm);
		int nprocs, myrank;
		MPI_Comm_rank(comm, &myrank);
		MPI_Comm_size(comm, &nprocs);
		mat.prepare_sequential_fill(5);
		
		for (size_t loc = 0; loc < mat.get_dim_local(); loc++)
		{
			size_t glob = loc + mat.first_row_on_node();
			if (glob == 0)
			{
				mat.sequential_fill(0, 4);
				mat.sequential_fill(1, -1);
				mat.sequential_fill(m, -1);
				// 4,-1 [...] -1
			}
			else if (myrank > 0 && myrank < m)
			{
				mat.sequential_fill(glob - 1, -1);
				mat.sequential_fill(glob, 4);
				mat.sequential_fill(glob + 1, -1);
				mat.sequential_fill(glob + m, -1);
				// -1, 4,-1 [...] - 1
			}
			else if (myrank >= m*(m - 1) && myrank != m*m - 1)
			{
				mat.sequential_fill(glob - m, -1);
				mat.sequential_fill(glob - 1, -1);
				mat.sequential_fill(glob, 4);
				mat.sequential_fill(glob + 1, -1);
				// -1 [...] -1 4 -1 
			}
			else if (myrank == m*m - 1)
			{
				mat.sequential_fill(glob - m, -1);
				mat.sequential_fill(glob - 1, -1);
				mat.sequential_fill(glob, 4);
				// -1 [...] -1 4
			}
			else
			{
				mat.sequential_fill(glob - m, -1);
				mat.sequential_fill(glob - 1, -1);
				mat.sequential_fill(glob, 4);
				mat.sequential_fill(glob + 1, -1);
				mat.sequential_fill(glob + m, -1);
				// -1 [...] -1 4 -1 [...] -1
			}
			mat.end_of_line();
		}
		return mat;
	}

	template <class Matrix, class Vector>
	class SolverBench
	{
		unsigned _node_min, _node_max, _m_min, _m_max;
		std::vector<std::vector<long> > _exec_times; // _exec_times[#nodes_idx][#m_idx]

	public:
		SolverBench(unsigned node_min, unsigned node_max, unsigned m_min, unsigned m_max) :
			_node_min(node_min),
			_node_max(node_max),
			_m_min(m_min),
			_m_max(m_max),
			_exec_times(node_max - node_min + 1, std::vector<long>(m_max - m_min + 1))
		{ }

		// Starte den Benchmark
		void run()
		{
			for (unsigned nodes = _node_min; nodes < _node_max; nodes++)
			{
				// sind wir an diesem benchmark beteiligt?
				int myglobalrank;
				MPI_Comm_rank(MPI_COMM_WORLD, &myglobalrank);
				if (myglobalrank >= nodes) continue;

				// konstruiere prozessgruppe mit #procs = nodes
				MPI_Barrier(MPI_COMM_WORLD);
				MPI_Comm pcomm;
				MPI_Group worldgroup, pgroup;
				MPI_Comm_group(MPI_COMM_WORLD, &worldgroup);
				int * ranks = new int[nodes];
				for (int i = 0; i < nodes; i++) ranks[i] = i;
				MPI_Group_incl(worldgroup, nodes, ranks, &pgroup);
				delete ranks;
				MPI_Comm_create_group(MPI_COMM_WORLD, pgroup, 0, &pcomm);
				int myrank = -1;
				MPI_Comm_rank(pcomm, &myrank);

				// messschleife
				for (unsigned m = _m_min; m < _m_max; m++)
				{
					// konstruiere matrix, startvektor und rechte seite
					Matrix mat = construct_model_matrix(m, pcomm);
					Vector rhs(m*m);
					rhs.fill_const(1.0);
					Vector res(m*m);
					res.clear();

					// Löser
					BiCgStabSolver<Matrix> solver(mat, rhs);

					// starte zeitmessung
					std::chrono::high_resolution_clock::time_point start;
					start = std::chrono::high_resolution_clock::now();

					// löser ausführen
					solver.solve();
					
					// speichere ergebnis
					_exec_times[nodes - _node_min][m - _m_min] = 
						std::chrono::duration_cast<std::chrono::milliseconds>(
						std::chrono::high_resolution_clock::now() - start).count();
				}

				// communicator und gruppe freigeben
				MPI_Comm_free(pcomm);
				MPI_Group_free(pgroup);
			}
		}

		// Drucke performace vs. #procs und performance/proc vs. #procs
		void print_results(std::ostream& out) const
		{
			out << "STRONG SCALING RESULTS:" << std::endl;
			for (unsigned m = _m_min; m < _m_max; m++)
				out << m*m << "\t";
			out << std::endl;
			for (unsigned nodes = _nodes_min; nodes < _node_max; nodes++)
			{
				// Drucke Zeile für #nodes [1. mwert 2.mwert ...]
				for (unsigned m = _m_min; m < _m_max; m++)
					out << 1.0/_exec_times[nodes - _modes_min][m - _m_min] << '\t';
				out << std::endl;
			}

			out << std::endl << "WEAK SCALING RESULTS:" << std::endl;
			for (unsigned m = _m_min; m < _m_max; m++)
				out << m*m << "\t";
			out << std::endl;
			for (unsigned nodes = _nodes_min; nodes < _node_max; nodes++)
			{
				// Drucke Zeile für #nodes [1. mwert 2.mwert ...]
				for (unsigned m = _m_min; m < _m_max; m++)
					out << 1.0/(nodes * _exec_times[nodes - _modes_min][m - _m_min]) << '\t';
				out << std::endl;
			}
		}
	};
}

#endif // __SOLVERBENCH_HPP_