#define NODEBUG
#define LOGGING_LEVEL 2

#include "mpihandler.hpp"
#include "slicedvector.hpp"
#include "distellpackmatrix.hpp"

#include <chrono>
#include <random>
#include <typeinfo>

int main(int nargs, char** args)
{
    using namespace std::chrono;
 
    constexpr size_t size = 1000000L;

    if(nargs != 3)
    {
        std::cout << "Usage: matvec [#nnz/row] [#samples]" << std::endl;
        exit(-1);
    }
    unsigned Ncol = atoi(args[1]);
    unsigned Nmean = atoi(args[2]);

    Icarus::DistEllpackMatrix<double> mat(size);
    mat.prepare_sequential_fill(Ncol);
    const size_t fron = mat.first_row_on_node();
    const size_t lron = fron + mat.get_dim_local();
    std::mt19937 eng;
    std::uniform_int_distribution<size_t> cdist(0,size-1);
    std::uniform_real_distribution<double> vdist(-1,1);
    for(size_t i = fron; i < lron; i++)
    {
        for(unsigned i=0; i< Ncol; i++)
            mat.sequential_fill(cdist(eng),vdist(eng));
        mat.end_of_row();
    }
    Icarus::SlicedVector<double> vec1(size), vec2(size);
    vec1.fill_const(vdist(eng));
    vec2.fill_const(vdist(eng));

    high_resolution_clock::time_point start =
        high_resolution_clock::now();

    for(unsigned i=0; i < Nmean; i++)
    {
	mat.mult_vec(vec1, vec2);
    }

    long long total_time = duration_cast<milliseconds>(
        high_resolution_clock::now() - start).count();

    double single_time = (double) total_time / Nmean;

    int myrank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Barrier(MPI_COMM_WORLD);
    if(myrank == 0)
    {
        std::cout << "Benchmark with " << nprocs << " processes finished." << std::endl;
        std::cout << "Total execution (wall) time: " << total_time << " ms" << std::endl
          << "Time per operation (1/" << Nmean << "): " << single_time << " ms" << std::endl; 
        std::cout << "ARM-CPU, random matrix with " << Ncol << " nnz/row." << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    return 0;
}
