#include <mpi.h>
#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "include/vmvr_help.hpp"
#include "include/roofline_help.hpp"
#include "include/timer.hpp"
using namespace std;

//------------------------------------------------------------------------------------------------/
//                                   APPLICATION SETTINGS
//------------------------------------------------------------------------------------------------/
//16384 - 65536 - 262144 - 1048576 - 4194304
#define dimension 16384
#define iteration 10000
bool get_overall = false;

template<typename type>
float invoke_gpu_time_dotproduct(type *vecx, type *vecy, type *placehold, int dim, int runs);

template<typename type>
float invoke_gpu_time_l2norm(type *vec, type *placehold, int dim, int runs);

template<typename type>
float invoke_gpu_time_reduce(type *placehold, type *result, int dim, int runs);

template<typename type>
void invoke_gpu_overall(type *vecx, type *vecy, type *result, int dim);

template<typename type>
void allocation(type **data, size_t size);

template <typename type>
void cleanup(type *data);


int main(int argc, char* argv[])
{
    int numblock = ceil((double)dim / 1024);
   
    double *vectorx_host = new double[dimension];
    double *vectory_host = new double[dimension];
    set_data(vectorx_host, dimension);
    set_data(vectory_host, dimension);

    Timer timer_overall;
    float elapsed_overall = 0.0;

//================================================================================================/
//									THE MAGIC HAPPENS HERE
//================================================================================================/
//------------------------------------------------------------------------------------------------/
//                                   Zeitmessung Overall
//------------------------------------------------------------------------------------------------/
    if (get_overall)
    {
        timer_overall.start();
        for (int r = 0; r < iteration; r++)
        {
            double *vectorx_dev = NULL;
            double *vectory_dev = NULL;
            double *placehold = NULL;

            allocation(&vectorx_dev, dimension);
            allocation(&vectory_dev, dimension);
            allocation(&placehold, numblock);

            copy_data(vectorx_host, vectorx_dev, dimension);
            copy_data(vectory_host, vectory_dev, dimension);

            invoke_gpu_overall(vectorx_dev, vectory_dev, placehold, dimension);

            cleanup(vectorx_dev);
            cleanup(vectory_dev);
            cleanup(placehold);
        }
        elapsed_overall = (timer_overall.stop()*1.0e3) / (float)iteration;
    }

    //------------------------------------------------------------------------------------------------/
    //                                   Skalarprodukt Block Result
    //------------------------------------------------------------------------------------------------/

    double *vectorx_dev = NULL;
    double *vectory_dev = NULL;
    double *placehold = NULL;
        allocation(&vectorx_dev, dimension);
        allocation(&vectory_dev, dimension);
        allocation(&placehold, numblock);
    copy_data(vectorx_host, vectorx_dev, dimension);
    copy_data(vectory_host, vectory_dev, dimension);
        float elapsed_kernel = invoke_gpu_time_dotproduct(vectorx_dev, vectory_dev, placehold, dimension, iteration);
    dotproduct_check_result(placehold, vectorx_host, vectory_host, dimension);
        cleanup(vectorx_dev);
        cleanup(vectory_dev);
        cleanup(placehold);

    double schalter = 0.0;
    //performance_dotproduct(dimension, elapsed_overall, elapsed_kernel, schalter, iteration);
  
    //------------------------------------------------------------------------------------------------/
    //                                   l2Norm Block Result
    //------------------------------------------------------------------------------------------------/

    double *vector_dev = NULL;
    double *placehold = NULL;
        allocation(&vector_dev, dimension);
        allocation(&placehold, numblock);
    copy_data(vectorx_host, vector_dev, dimension);
        elapsed_kernel = invoke_gpu_time_l2norm(vector_dev, placehold, dimension, iteration);
    l2norm_check_result_(placehold, vector_host, dimension);

    //performance_l2norm(dimension, elapsed_overall, elapsed_kernel, schalter, iteration);
    
    //COMPUTE TIME OF REDUCE KERNEL
    double *result = NULL;
    allocation(&result, 1);
    elapsed_kernel = invoke_gpu_time_reduce(placehold,reduce,numblocks);
    reduce_check_result(result, placehold, numblock);
    
    //performance_l2norm(dimension, elapsed_overall, elapsed_kernel, schalter, iteration);
    
    cleanup(vector_dev);
    cleanup(placehold);
    cleanup(result);

    //================================================================================================/
    //                                         Evaluieren
    //================================================================================================/

    delete[] vectorx_host;
    delete[] vectory_host;

    return 0;
}

