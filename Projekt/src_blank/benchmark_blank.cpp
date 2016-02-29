#include <mpi.h>
#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "include/benchmark_help.hpp"
#include "include/timer.hpp"
using namespace std;
enum methods { unified_vs_zero, kernel_vs_kernel, kernel_vs_cpu };  //choose a method
enum version { kernel_standart, kernel_shared, kernel_advanced };   //keep your kernels in the same order as they are in the switch in your gpu_xx_call
enum memory_opt { unified, zero };                                  //choose a method of memory usage(for method k_vs_k and k_vs_cpu)
//================================================================================================/
//									GLOBAL SETTINGS!
//================================================================================================/
int method = kernel_vs_kernel;               
int version_first = kernel_standart;        
//Für Kernel vs Kernel
int version_second = kernel_standart;      
//Für Kernel vs Kernel und Kernel vs CPU
int memory_option = zero;

//------------------------------------------------------------------------------------------------/
//                                   APPLICATION SETTINGS
//------------------------------------------------------------------------------------------------/


void print_p();

template <typename type>
void cpu_ax(type *pointer);

template<typename type>
void performance(float time_ku, float time_ou, float time_kz, float time_oz, int runs, type schalter, int meth, int ver_first, int ver_second, int mem_option);

template<typename Scalar>
void gpu_ax_overall(Scalar *pointer, int version, int mem_option);

template<typename Scalar>
float gpu_ax_time(Scalar *pointer, int runs, int version, int mem_option);

template<typename Scalar>
void allocation(Scalar **pointer, int mem_option);

template <typename Scalar>
void cleanup(Scalar *pointer, int method);

int main(int argc, char* argv[])
{
    //Array zur Zeitmessung
    //Generiere data/Indices Int-Array sowie fvec Int Array
    

    Timer timer_overall,timer_cpu;

//================================================================================================/
//									THE MAGIC HAPPENS HERE
//================================================================================================/
//------------------------------------------------------------------------------------------------/
//                                Zeitmessung Overall Teil 1
//------------------------------------------------------------------------------------------------/
    if (method == unified_vs_zero) { memory_option = unified; }

    timer_overall.start();
    for(int r = 0;r<iteration;r++)
    {
        float *pointer_first = NULL;
      
        allocation(&pointer_first, memory_option);
        gpu_ax_overall(pointer_first, version_first, memory_option);
        cleanup(pointer_first, memory_option);
    }
    float elapsed_first_overall = timer_overall.stop() / (float)iteration;

//------------------------------------------------------------------------------------------------/
//                                Zeitmessung Overall Teil 2
//------------------------------------------------------------------------------------------------/
    if (method == unified_vs_zero) { memory_option = zero; version_second = version_first; }

    timer_overall.start();
    for (int r = 0; r<iteration; r++)
    {
        float *pointer_second = NULL;

        if (method != kernel_vs_cpu)
        {
            allocation(&pointer_second, memory_option);
            gpu_ax_overall(pointer_second, version_second, memory_option);
            cleanup(pointer_second, memory_option);
        }
        else//CPU Zeitmessung
        {
            
            cpu_ax(pointer_second);
            cleanup(pointer_second, 2);
        }
    }
    float elapsed_second_overall = timer_overall.stop() / (float)iteration;


//------------------------------------------------------------------------------------------------/
//                                Zeitmessung Kernel Teil 1                   
//------------------------------------------------------------------------------------------------/
    if (method == unified_vs_zero) { memory_option = unified; }
    
    float *pointer_first = NULL;

    allocation(&pointer_first, memory_option);

    //=========================================//Hier muss vielleicht die Zeitmessung innerhalb der aufgerufenen Funktion stattfinden
    float elapsed_first_kernel = 
        gpu_ax_time(pointer_first, iteration, version_first, memory_option);
    //=========================================//

    cleanup(pointer_first, memory_option);
 
 //------------------------------------------------------------------------------------------------/
 //                                Zeitmessung Kernel Teil 2                   
 //------------------------------------------------------------------------------------------------/
    if (method == unified_vs_zero) { memory_option = zero; version_second = version_first; }
    
    float *pointer_second = NULL;
    float elapsed_second_kernel = 0.0;

    if (method != kernel_vs_cpu)
    {
        allocation(&pointer_second, memory_option);
        
        //=========================================//Hier muss vielleicht die Zeitmessung innerhalb der aufgerufenen Funktion stattfinden
        elapsed_second_kernel =
            gpu_ax_time(pointer_second, iteration, version_second, memory_option);
        //=========================================//
        
        cleanup(pointer_second, memory_option);
    }
    else//CPU Zeitmessung
    {
        
        
        //=========================================//
        timer_cpu.start();
        for (int r = 0; r < iteration; r++)
        {
            cpu_ax(pointer_second);
        }
        elapsed_second_kernel = timer_cpu.stop()*1.0e3;
        //=========================================//

        cleanup(pointer_second, 2);
    }
    
   
//================================================================================================/
//                                         Evaluieren
//================================================================================================/

    float schalter = 0.0;
    performance(elapsed_first_kernel, elapsed_first_overall, elapsed_second_kernel, elapsed_second_overall, iteration, schalter, method, version_first, version_second, memory_option);
  
    delete[] pointer_host;

    return 0;
}

template <typename type>
void cpu_ax(type *pointer)
{

}
template void cpu_ax<int>(int *pointer);
template void cpu_ax<float>(float *pointer);
template void cpu_ax<double>(double *pointer);