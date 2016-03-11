#include <mpi.h>
#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "include/copy_help.hpp"
#include "include/roofline_help.hpp"
#include "include/timer.hpp"
using namespace std;

//------------------------------------------------------------------------------------------------/
//                                   APPLICATION SETTINGS
//------------------------------------------------------------------------------------------------/
#define dimension 32768
#define iteration 1000

template<typename type>
float invoke_gpu_time(type *vecin, type scalar, type *vecout, int dim, int runs);

template<typename type>
void invoke_gpu_overall(type *vecin, type scalar, type *vecout, int dim);

template<typename type>
void allocation(type **data, size_t size);

template <typename type>
void cleanup(type *data);


int main(int argc, char* argv[])
{
   
    double *vec_host = new double[dimension];
    double *scalar_host = new double[1];
    set_data(vec_host, dimension);
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
        double *vecin_dev = NULL;
        double *vecout_dev = NULL;

        allocation(&vecin_dev, dimension);
        allocation(&vecout_dev, dimension);
        
        copy_data(vec_host, vecin_dev, dimension);

        invoke_gpu_overall(vecin_dev, scalar_host[0], vecout_dev, dimension);

        cleanup(vecin_dev);
        cleanup(vecout_dev);
    }
    float elapsed_overall = timer_overall.stop() / (float)iteration;


//------------------------------------------------------------------------------------------------/
//                                   Zeitmessung Kernel
//------------------------------------------------------------------------------------------------/
    
    double *vecin_dev = NULL;
    double *vecout_dev = NULL;

    allocation(&vecin_dev, dimension);
    allocation(&vecout_dev, dimension);

    copy_data(vec_host, vecin_dev, dimension);

    //=========================================//Hier muss vielleicht die Zeitmessung innerhalb der aufgerufenen Funktion stattfinden
    float elapsed_kernel = invoke_gpu_time(vecin_dev, scalar_host[0], vecout_dev, dimension, iteration);
        //>>>KERNEL<<<
    //=========================================//

    copy_check_result(vec_host, scalar_host[0], vecout_dev, dimension);
    cleanup(vecin_dev);
    cleanup(vecout_dev);

 
//================================================================================================/
//                                         Evaluieren
//================================================================================================/

    double schalter = 0.0;
    performance(dimension, elapsed_overall, elapsed_kernel, schalter);
  
    delete[] vec_host;
    delete[] scalar_host;

    return 0;
}

