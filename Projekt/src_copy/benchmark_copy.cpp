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
int method = unified_vs_zero;
int version_first = kernel_standart;        
//Für Kernel vs Kernel
int version_second = kernel_standart;      
//Für Kernel vs Kernel und Kernel vs CPU
int memory_option = zero;

//------------------------------------------------------------------------------------------------/
//                                   APPLICATION SETTINGS
//------------------------------------------------------------------------------------------------/
#define dimlocal 16384
#define iteration 1000


void print_p();

//template <typename type>
//void cpu_ax(type *pointer);

template<typename type>
void performance(float time_ku, float time_ou, float time_kz, float time_oz, int runs, type schalter, int meth, int ver_first, int ver_second, int mem_option, int dim_local);

template<typename Scalar>
void gpu_copy_overall(Scalar *one, Scalar *two, int dim_local, int version, int mem_option);

template<typename Scalar>
float gpu_copy_time(Scalar *one, Scalar *two, int dim_local, int runs, int version, int mem_option);

/*DONE*/template<typename Scalar>
void allocation(Scalar **vecin, Scalar **vectwo, int dim_local, int mem_option);

/*DONE*/template <typename Scalar>
void cleanup(Scalar *one, Scalar *two, int method);

int main(int argc, char* argv[])
{
    //Array zur Zeitmessung
    //Generiere data/Indices Int-Array sowie fvec Int Array
    float *vecin_host = new float[dimlocal];
    vec_float_one(vector_host, dimlocal);

    Timer timer_overall,timer_cpu;

//================================================================================================/
//									THE MAGIC HAPPENS HERE
//================================================================================================/
//------------------------------------------------------------------------------------------------/
//                                Zeitmessung Overall Teil 1
//------------------------------------------------------------------------------------------------/
    if (method == unified_vs_zero) 
    { 
        memory_option = unified; 
    }

    timer_overall.start();
    for(int r = 0;r<iteration;r++)
    {
        float *vecin_first = NULL;
        float *vecout_first = NULL;
      
        allocation(&vecin_first, &vecout_first, dimlocal, memory_option);
        set_values(vecin_host, vecin_first, dimlocal);
        
        gpu_dotproduct_overall(vecin_first, vecout_first, dimlocal, version_first, memory_option);
        cleanup(vecin_first, vecout_first, memory_option);
    }
    float elapsed_first_overall = timer_overall.stop() / (float)iteration;


//------------------------------------------------------------------------------------------------/
//                                Zeitmessung Overall Teil 2
//------------------------------------------------------------------------------------------------/
    if (method == unified_vs_zero) 
    { 
        memory_option = zero; 
        version_second = version_first; 
    }
    timer_overall.start();
    for (int r = 0; r<iteration; r++)
    {
        float *vecin_second = NULL;
        float *vecout_second = NULL;

        if (method != kernel_vs_cpu)
        {
            allocation(&vecin_second, &vecout_second, dimlocal, memory_option);
            set_values(vecin_host, vecin_second, dimlocal);
            gpu_dotproduct_overall(vecin_second, vecout_second, dimlocal, version_second, memory_option);
            cleanup(vecin_second, vecout_second, memory_option);
        }
        else//CPU Zeitmessung
        {
            
            //cpu_ax(pointer_second);
            //cleanup(pointer_second, 2);
        }
    }
    float elapsed_second_overall = timer_overall.stop() / (float)iteration;


//------------------------------------------------------------------------------------------------/
//                                Zeitmessung Kernel Teil 1                   
//------------------------------------------------------------------------------------------------/
    if (method == unified_vs_zero) { memory_option = unified; }
    float *vecin_first = NULL;
    float *vecout_first = NULL;

    allocation(&vecin_first, &vecout_first, dimlocal, memory_option);
    set_values(vecin_host, vecin_first, dimlocal);

    //=========================================//Hier muss vielleicht die Zeitmessung innerhalb der aufgerufenen Funktion stattfinden
    float elapsed_first_kernel = 
        gpu_dotproduct_time(vecin_first, vecout_first, dimlocal, iteration, version_first, memory_option);
    //=========================================//

    check_vecout_dotproduct(vecout_first, vecin_host, dimlocal, '1');
    cleanup(vecin_first, vecout_first, memory_option);
 
 //------------------------------------------------------------------------------------------------/
 //                                Zeitmessung Kernel Teil 2                   
 //------------------------------------------------------------------------------------------------/
    if (method == unified_vs_zero) { memory_option = zero; version_second = version_first; }
    float *vecin_second = NULL;
    float *vecout_second = NULL;

    float elapsed_second_kernel = 0.0;

    if (method != kernel_vs_cpu)
    {
        allocation(&vecin_second, &vecout_second, dimlocal, memory_option);
        set_values(vecin_host, vecin_second, dimlocal);
        
        //=========================================//Hier muss vielleicht die Zeitmessung innerhalb der aufgerufenen Funktion stattfinden
        elapsed_second_kernel =
            gpu_dotproduct_time(vecin_second, vecout_second, dimlocal, iteration, version_second, memory_option);
        //=========================================//
        
        check_vecout_dotproduct(vecout_second, vecin_host, dimlocal, '2');
        cleanup(vecin_second, vecout_second, memory_option);
    }
    else//CPU Zeitmessung
    {
        
        
        //=========================================//
        //timer_cpu.start();
        //for (int r = 0; r < iteration; r++)
        //{
            //cpu_ax(pointer_second);
       // }
        //elapsed_second_kernel = timer_cpu.stop()*1.0e3;
        //=========================================//

        //cleanup(pointer_second, 2);
    }
    
 
//================================================================================================/
//                                         Evaluieren
//================================================================================================/

    float schalter = 0.0;
    performance(elapsed_first_kernel, elapsed_first_overall, elapsed_second_kernel, elapsed_second_overall, iteration, schalter, method, version_first, version_second, memory_option, dimlocal);
  
    delete[] vecin_host;
    delete[] vecout_host;

    return 0;
}

/*template <typename type>
void cpu_ax(type *pointer)
{

}
template void cpu_ax<int>(int *pointer);
template void cpu_ax<float>(float *pointer);
template void cpu_ax<double>(double *pointer);*/