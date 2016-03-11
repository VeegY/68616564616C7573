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
float invoke_gpu_time(int runs);

template<typename type>
void invoke_gpu_overall();

template<typename type>
void allocation(type **data, size_t size);

template <typename type>
void cleanup(type *data);


int main(int argc, char* argv[])
{
   
    double *data_host = new double[dimension];
    set_data(data_host, dimension);

    Timer timer_overall,timer_cpu;

//================================================================================================/
//									THE MAGIC HAPPENS HERE
//================================================================================================/
//------------------------------------------------------------------------------------------------/
//                                   Zeitmessung Overall
//------------------------------------------------------------------------------------------------/

    timer_overall.start();
    for(int r = 0;r<iteration;r++)
    {
        double *data_dev = NULL;
        allocation(&data_dev, dimension);
        copy_data(data_host, data_dev, dimension);

        >>>KERNEL<<<<

        cleanup(data_dev);
    }
    float elapsed_overall = timer_overall.stop() / (float)iteration;


//------------------------------------------------------------------------------------------------/
//                                   Zeitmessung Kernel
//------------------------------------------------------------------------------------------------/
    
    double *data_dev = NULL;
    allocation(&data_dev, dimension);
    copy_data(data_host, data_dev, dimension);

    //=========================================//Hier muss vielleicht die Zeitmessung innerhalb der aufgerufenen Funktion stattfinden
    float elapsed_kernel =
        >>>KERNEL<<<
    //=========================================//

    //check_result!
    cleanup(data_dev);

 
//================================================================================================/
//                                         Evaluieren
//================================================================================================/

    double schalter = 0.0;
    //performance(elapsed_kernel, elapsed_overall, schalter);
  
    delete[] data_host;

    return 0;
}

