/*
* Autor       : David Schneider
* Datum       : 9.4.2016
* Beschreibung: Demo zur Skalierbarkeit der MV-Multiplikation
*               von Bandmatrizen mit NN-Kommunikation
*/

#define NDEBUG

#include <chrono>
#include <iostream>
#include <omp.h>
#include <cuda_runtime.h>

#include "bandscal.h"

void print_usage();

void parse_args(int nargs, char** args, arch_t& arch, bool& dp, int& m);

int main(int nargs, char** args)
{
    // cmd line
    arch_t arch;
    bool dp;
    int m;
    omp_set_dynamic(0);
    parse_args(nargs, args, arch, dp, m);

    const int nsamples = 100;

    // init
    MPI_Init(&nargs, &args);
    int myrank, nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    cublasHandle_t cublas_handle;
    cublasStatus_t err1 = cublasCreate(&cublas_handle);
    cusparseHandle_t cusp_handle;
    cusparseStatus_t err2 = cusparseCreate(&cusp_handle);
    if(err1 != CUBLAS_STATUS_SUCCESS || err2 != CUSPARSE_STATUS_SUCCESS)
    {
      std::cerr << "Failed to initialize cuBLAS / cuSPARSE." << std::endl;
      std::exit(-1);
    }
    int nthreads;
#pragma omp parallel
    {
        nthreads = omp_get_num_threads();
    }
    for (int i = 0; i < nprocs; i++)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        if (myrank != i) continue;
        std::cout << "Process " << myrank << " of " << nprocs 
            << " has " << nthreads << " thread(s) and "
            << " successfully initialized MPI, cuBLAS and cuSPARSE." << std::endl;
    }
    
    // Benchmark
    if(myrank == 0) std::cout << "Setting up benchmark..." << std::endl;
    try
    {
        if (dp) // double precision
        {
            BCsrMatrix<double> mat = construct_model_matrix<double>(m, arch, cusp_handle);
            BVector<double> x(m*m, m, arch, cublas_handle), y(m*m, m, arch, cublas_handle);
            x.fill_with(1.0);
            if(myrank == 0) std::cout << "Starting bechmark." << std::endl;
            std::chrono::high_resolution_clock::time_point start =
                std::chrono::high_resolution_clock::now();
            for(int i=0; i <nsamples; i++) mat.spmv(x,y);
            std::chrono::duration<double, std::milli> elapsed =
                std::chrono::high_resolution_clock::now() - start;
            if (myrank == 0)
                std::cout << "Elapsed time for matvec (ms):" << std::endl << elapsed.count()/nsamples << std::endl;
            y.print(std::cout);
        }
        else // single precision
        {
            BCsrMatrix<float> mat = construct_model_matrix<float>(m, arch, cusp_handle);
            BVector<float> x(m*m, m, arch, cublas_handle), y(m*m, m, arch, cublas_handle);
            x.fill_with(1.0);
            if(myrank == 0) std::cout << "Starting bechmark." << std::endl;
            std::chrono::high_resolution_clock::time_point start =
                std::chrono::high_resolution_clock::now();
            for(int i=0; i <nsamples; i++) mat.spmv(x,y);
            std::chrono::duration<double, std::milli> elapsed =
                std::chrono::high_resolution_clock::now() - start;
            if (myrank == 0)
                std::cout << "Elapsed time for matvec (ms):" << std::endl  << elapsed.count() << std::endl;
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "Process " << myrank << " caught an exception: " << e.what() << std::endl;
    }
    catch (...)
    {
        std::cerr << "Process " << myrank << " caught an unexpected exception." << std::endl;
    }

    // cleanup
    cusparseDestroy(cusp_handle);
    cublasDestroy(cublas_handle);
    MPI_Finalize();
    return 0;
}

void print_usage()
{
    std::cout << ">>> Usage: mvb [arch/threads] [prec] [size]" << std::endl
        << "arch: {'G','1','2','3','4'} for GPU or up to 4 threads on ARM" << std::endl
        << "prec: {'s','d'} for single or double precision" << std::endl
        << "size: specify parameter m in model matrix" << std::endl;
}

void parse_args(int nargs, char** args, arch_t& arch, bool& dp, int& m)
{
    if (nargs != 4)
    {
        print_usage();
        exit(-1);
    }
    switch (args[1][0])
    {
    case 'g':
    case 'G':
        arch = ARCH_GPU;
        omp_set_num_threads(1);
        break;
    case '1':
        arch = ARCH_CPU;
        omp_set_num_threads(1);
        break;
    case '2':
        arch = ARCH_CPU;
        omp_set_num_threads(2);
        break;
    case '3':
        arch = ARCH_CPU;
        omp_set_num_threads(3);
        break;
    case '4':
        arch = ARCH_CPU;
        omp_set_num_threads(4);
        break;
    default:
        print_usage();
        exit(-1);
    }

    switch (args[2][0])
    {
    case 'S':
    case 's':
        dp = false;
        break;
    case 'D':
    case 'd':
        dp = true;
        break;
    default:
        print_usage();
        exit(-1);
    }

    m = atoi(args[3]);
}
