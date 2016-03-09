#ifndef __SOLVERBENCH_HPP_
#define __SOLVERBENCH_HPP_

#include <iostream>
#include <vector>
#include <mpi.h>
#include <chrono>
#include "bicgstabsolver.hpp"
#include "distellpackmatrix.hpp"
#include "logger.hpp"
#include "utility.hpp"

namespace Icarus
{
	// ACHTUNG: Nur die Prozesse, die in comm enthalten sind, dürfen diese Funktion aufrufen.
	template <class DistMatrix>
	static DistMatrix construct_model_matrix(unsigned m, MPI_Comm comm = MPI_COMM_WORLD)
	{
		if (comm == MPI_COMM_NULL)
		{
			std::cerr << "construct_model_matrix von ungueltigem comm aus aufgerufen!" << std::endl;
			throw;
		}

		DistMatrix mat(m*m, comm);
		mat.prepare_sequential_fill(5);
		
		for (size_t loc = 0; loc < mat.get_dim_local(); loc++)
		{
			const size_t glob = loc + mat.first_row_on_node();
			if (glob == 0)
			{
				mat.sequential_fill(0, 4);
				mat.sequential_fill(1, -1);
				mat.sequential_fill(m, -1);
				// 4,-1 [...] -1
			}
			else if (glob > 0 && glob < m)
			{
				mat.sequential_fill(glob - 1, -1);
				mat.sequential_fill(glob, 4);
				mat.sequential_fill(glob + 1, -1);
				mat.sequential_fill(glob + m, -1);
				// -1, 4,-1 [...] - 1
			}
			else if (glob >= m*(m - 1) && glob != m*m - 1)
			{
				mat.sequential_fill(glob - m, -1);
				mat.sequential_fill(glob - 1, -1);
				mat.sequential_fill(glob, 4);
				mat.sequential_fill(glob + 1, -1);
				// -1 [...] -1 4 -1 
			}
			else if (glob == m*m - 1)
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
			mat.end_of_row();
		}
		return mat;
	}

	template <class Matrix>
	class SolverBench
	{
		std::vector<unsigned> _nlist, _mlist;
		std::vector<std::vector<long> > _exec_times; // _exec_times[#nodes_idx][#m_idx]

	public:
		SolverBench(unsigned node_min, unsigned node_max, unsigned m_min, unsigned m_max) :
			_nlist(node_max - node_min + 1),
			_mlist(m_max - m_min + 1),
			_exec_times(node_max - node_min + 1, std::vector<long>(m_max - m_min + 1))
		{
			int nnodes;
			MPI_Comm_size(MPI_COMM_WORLD, &nnodes);
			if(nnodes < (int)node_max)
				LOG_ERROR("Not enough nodes available for ths benchmark (", nnodes, " provided, ", node_max, " required).");
			// listen füllen
			for(unsigned n = node_min; n <= node_max; n++) _nlist[n-node_min] = n;	
			for(unsigned m = m_min; m <= m_max; m++) _mlist[m-m_min] = m;	
		}

		SolverBench(const std::vector<unsigned>& nlist, const std::vector<unsigned>& mlist) :
			_nlist(nlist),
			_mlist(mlist),
			_exec_times(nlist.size(), std::vector<long>(mlist.size()))
		{
			int nnodes;
			MPI_Comm_size(MPI_COMM_WORLD, &nnodes);
			for(unsigned n : nlist) if((int)n > nnodes)
				LOG_ERROR("Not enough nodes available for ths benchmark (", nnodes, " provided, >=", n, " required).");
		}

		// Starte den Benchmark
		void run()
		{
			for (unsigned nctr = 0; nctr < _nlist.size(); nctr++)
			{
				const unsigned nodes = _nlist[nctr];
				MPI_Barrier(MPI_COMM_WORLD);
				LOG_INFO("Now starting benchmark with ", nodes, " node(s).");

				// konstruiere prozessgruppe mit #procs = nodes
				MPI_Comm pcomm;
				MPI_Group worldgroup, pgroup;
				MPI_Comm_group(MPI_COMM_WORLD, &worldgroup);
				int * ranks = new int[nodes];
				for (int i = 0; i < (int)nodes; i++) ranks[i] = i;
				MPI_Group_incl(worldgroup, nodes, ranks, &pgroup);
				delete ranks;
				MPI_Comm_create_group(MPI_COMM_WORLD, pgroup, 0, &pcomm);
				int myrank = -1;
				MPI_Comm_rank(pcomm, &myrank);
				LOG_DEBUG("Process group successfully created.");

				// sind wir an diesem benchmark beteiligt?
				int myglobalrank;
				MPI_Comm_rank(MPI_COMM_WORLD, &myglobalrank);
				if (myglobalrank < (int)nodes)
				{
				
				// messschleife
				for (unsigned mctr = 0; mctr < _mlist.size(); mctr++)
				{
					const unsigned m = _mlist[mctr]; 
					// konstruiere matrix, startvektor und rechte seite
					Matrix mat = construct_model_matrix<Matrix>(m, pcomm);
					typename Matrix::VectorType rhs(m*m, pcomm);
					rhs.fill_const(1.0);
					typename Matrix::VectorType res(m*m, pcomm);
					res.clear();

					// Löser
					BiCgStabSolver<Matrix> solver(mat, rhs);

					// starte zeitmessung
					std::chrono::high_resolution_clock::time_point start;
					start = std::chrono::high_resolution_clock::now();

					// löser ausführen
					solver.solve(res);
					
					// speichere ergebnis
					_exec_times[nctr][mctr] = 
						std::chrono::duration_cast<std::chrono::milliseconds>(
						std::chrono::high_resolution_clock::now() - start).count();
					}
				}

				// communicator und gruppe freigeben
				MPI_Comm_free(&pcomm);
				MPI_Group_free(&pgroup);
			}
		}

		// Drucke performace vs. #procs und performance/proc vs. #procs
		void print_results(std::ostream& out, int rank = 0) const
		{
                        int myrank;
			MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
                        if(myrank != rank) return;

			out << "STRONG SCALING RESULTS:" << std::endl;
			for (unsigned mctr = 0; mctr < _mlist.size(); mctr++)
				out << _mlist[mctr]*_mlist[mctr] << "\t";
			out << std::endl;
			for (unsigned nctr = 0; nctr < _nlist.size(); nctr++)
			{
				// Drucke Zeile für #nodes [1. mwert 2.mwert ...]
				for (unsigned mctr = 0; mctr < _mlist.size(); mctr++)
					out << 1.0/_exec_times[nctr][mctr] << '\t';
				out << std::endl;
			}

			out << std::endl << "WEAK SCALING RESULTS:" << std::endl;
			for (unsigned mctr = 0; mctr < _mlist.size(); mctr++)
				out << _mlist[mctr]*_mlist[mctr] << "\t";
			out << std::endl;
			for (unsigned nctr = 0; nctr < _nlist.size(); nctr++)
			{
				// Drucke Zeile für #nodes [1. mwert 2.mwert ...]
				for (unsigned mctr = 0; mctr < _mlist.size(); mctr++)
					out << 1.0/(_nlist[nctr] * _exec_times[nctr][mctr]) << '\t';
				out << std::endl;
			}
		}
	};
}

#endif // __SOLVERBENCH_HPP_
