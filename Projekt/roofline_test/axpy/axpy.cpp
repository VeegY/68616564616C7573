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
#define dimension 1024
#define iteration 1000

template<typename type>
float invoke_gpu_time(type scalar, type *vector_x, type *vector_y, type *result, int dim, int runs);

template<typename type>
void invoke_gpu_overall(type scalar, type *vector_x, type *vector_y, type *result, int dim);

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
    set_data(scalar, 1);

    Timer timer_overall,timer_cpu;

//================================================================================================/
//									THE MAGIC HAPPENS HERE
//================================================================================================/
//------------------------------------------------------------------------------------------------/
//                                   Zeitmessung Overall
//------------------------------------------------------------------------------------------------/
    if (method == unified_vs_zero) 
    { 
        memory_option = unified; 
    }

    timer_overall.start();
    for(int r = 0;r<iteration;r++)
    {
        double *scalar = NULL;
        double *vectorx = NULL;
        double *vectory = NULL;
        double *result = NULL;
      
        allocation(&scalar, 1);
        allocation(&vectorx, dimension);
        allocation(&vectory, dimension);
        allocation(&result, dimension);

        copy_data(scalar_host, scalar, 1);
        copy_data(vectorx_host, vectorx, dimension);
        copy_data(vectory_host, vectory, dimension);
        
        >>>KERNEL<<<<

        cleanup(scalar);
        cleanup(vectorx);
        cleanup(vectory);
    }
    float elapsed_overall = timer_overall.stop() / (float)iteration;


//------------------------------------------------------------------------------------------------/
//                                   Zeitmessung Kernel
//------------------------------------------------------------------------------------------------/
    
    double *scalar = NULL;
    double *vectorx = NULL;
    double *vectory = NULL;
    double *result = NULL;

    allocation(&scalar, 1);
    allocation(&vectorx, dimension);
    allocation(&vectory, dimension);
    allocation(&result, dimension);

    copy_data(scalar_host, scalar, 1);
    copy_data(vectorx_host, vectorx, dimension);
    copy_data(vectory_host, vectory, dimension);

    //=========================================//Hier muss vielleicht die Zeitmessung innerhalb der aufgerufenen Funktion stattfinden
    float elapsed_kernel =
        >>>KERNEL<<<
    //=========================================//


 
//================================================================================================/
//                                         Evaluieren
//================================================================================================/

    double schalter = 0.0;
    performance(elapsed_kernel, elapsed_overall, dimension, schalter);
  
    delete[] vectorx_host;
    delete[] vectory_host;
    delete[] scalar;

    return 0;
}

