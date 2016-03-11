#include <mpi.h>
#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "include/axpy_help.hpp"
#include "include/roofline_help.hpp"
#include "include/timer.hpp"
using namespace std;

//------------------------------------------------------------------------------------------------/
//                                   APPLICATION SETTINGS
//------------------------------------------------------------------------------------------------/
#define dimension 32768
#define iteration 1000

template<typename type>
float invoke_gpu_time(type *vector_x, type scalar, type *vector_y, type *result, int dim, int runs);

template<typename type>
void invoke_gpu_overall(type *vector_x, type scalar, type *vector_y, type *result, int dim);

template<typename type>
void allocation(type **data, size_t size);

template <typename type>
void cleanup(type *data);

int main(int argc, char* argv[])
{
   
    double *vectorx_host = new double[dimension];
    double *vectory_host = new double[dimension];
    double *scalar_host = new double[1];

    set_data(vectorx_host, dimension);
    set_data(vectory_host, dimension);
    set_data(scalar_host, 1);

    Timer timer_overall;

//================================================================================================/
//									THE MAGIC HAPPENS HERE
//================================================================================================/
//------------------------------------------------------------------------------------------------/
//                                   Zeitmessung Overall
//------------------------------------------------------------------------------------------------/

    timer_overall.start();
    for(int r = 0;r<iteration;r++)
    {

        double *vectorx = NULL;
        double *vectory = NULL;
        double *result = NULL;
      
        allocation(&vectorx, dimension);
        allocation(&vectory, dimension);
        allocation(&result, dimension);

        copy_data(vectorx_host, vectorx, dimension);
        copy_data(vectory_host, vectory, dimension);
        
        invoke_gpu_overall(vectorx, scalar_host[0], vectory, result, dimension);

        cleanup(vectorx);
        cleanup(vectory);
    }
    float elapsed_overall = (timer_overall.stop()*1.0e3) / (float)iteration;


//------------------------------------------------------------------------------------------------/
//                                   Zeitmessung Kernel
//------------------------------------------------------------------------------------------------/
    
    double *vectorx = NULL;
    double *vectory = NULL;
    double *result = NULL;

    allocation(&vectorx, dimension);
    allocation(&vectory, dimension);
    allocation(&result, dimension);

    copy_data(vectorx_host, vectorx, dimension);
    copy_data(vectory_host, vectory, dimension);

    //=========================================//Hier muss vielleicht die Zeitmessung innerhalb der aufgerufenen Funktion stattfinden
    float elapsed_kernel = invoke_gpu_time(vectorx, scalar_host[0], vectory, result, dimension, iteration);
    //>>>KERNEL<<<
    //=========================================//

    axpy_check_result_(result, vectorx_host, scalar_host[0], vectory_host, dimension);
    cleanup(vectorx);
    cleanup(vectory);

 
//================================================================================================/
//                                         Evaluieren
//================================================================================================/

    double schalter = 0.0;
    performance(dimension, elapsed_overall, elapsed_kernel, schalter, iteration);
  
    delete[] vectorx_host;
    delete[] vectory_host;
    delete[] scalar_host;

    return 0;
}

