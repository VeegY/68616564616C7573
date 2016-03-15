#include <mpi.h>
#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "include/dotproduct_help.hpp"
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
float invoke_gpu_time(type *vecx, type *vecy, type *result, int dim, int runs);

template<typename type>
void invoke_gpu_overall(type *vecx, type *vecy, type *result, int dim);

template<typename type>
void allocation(type **data, size_t size);

template <typename type>
void cleanup(type *data);


int main(int argc, char* argv[])
{
   
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
            double *result = NULL;

            allocation(&vectorx_dev, dimension);
            allocation(&vectory_dev, dimension);
            allocation(&result, 1);

            copy_data(vectorx_host, vectorx_dev, dimension);
            copy_data(vectory_host, vectory_dev, dimension);

            invoke_gpu_overall(vectorx_dev, vectory_dev, result, dimension);

            cleanup(vectorx_dev);
            cleanup(vectory_dev);
            cleanup(result);
        }
        elapsed_overall = (timer_overall.stop()*1.0e3) / (float)iteration;
    }

//------------------------------------------------------------------------------------------------/
//                                   Zeitmessung Kernel
//------------------------------------------------------------------------------------------------/
    
    double *vectorx_dev = NULL;
    double *vectory_dev = NULL;
    double *result = NULL;

    allocation(&vectorx_dev, dimension);
    allocation(&vectory_dev, dimension);
    allocation(&result, 1);

    copy_data(vectorx_host, vectorx_dev, dimension);
    copy_data(vectory_host, vectory_dev, dimension);

    //=========================================//Hier muss vielleicht die Zeitmessung innerhalb der aufgerufenen Funktion stattfinden
    float elapsed_kernel = invoke_gpu_time(vectorx_dev, vectory_dev, result, dimension, iteration);
    //>>>KERNEL<<<
    //=========================================//

    
    dotproduct_check_result_(result, vectorx_host, vectory_host, dimension);
    cleanup(vectorx_dev);
    cleanup(vectory_dev);
    cleanup(result);
 
//================================================================================================/
//                                         Evaluieren
//================================================================================================/

    double schalter = 0.0;
    performance(dimension, elapsed_overall, elapsed_kernel, schalter, iteration);
  
    delete[] vectorx_host;
    delete[] vectory_host;

    return 0;
}

