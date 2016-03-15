#include <mpi.h>
#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "include/matvec_help.hpp"
#include "include/roofline_help.hpp"
#include "include/timer.hpp"
using namespace std;

//                                   APPLICATION SETTINGS
//------------------------------------------------------------------------------------------------/
//16384 - 65536 - 262144 - 1048576 - 4194304
#define dimension 65536
#define maxrowlength 7
#define iteration 10000
bool get_overall = false;
//------------------------------------------------------------------------------------------------/

template<typename type>
float invoke_gpu_time(type *vector, type *data, int * indices, type *result, int dim, int max_row_length, int runs);

template<typename type>
void invoke_gpu_overall(type *vector, type *data, int * indices, type *result, int dim, int max_row_length);

template<typename type>
void allocation(type **data, size_t size);

template <typename type>
void cleanup(type *data);


int main(int argc, char* argv[])
{
   
    double *vector_host = new double[dimension];
    double *data_host = new double[dimension*maxrowlength];
    int *indices_host = new int[dimension*maxrowlength];
    
    set_data(vector_host, dimension);
    ellpack_fill_seven_diagonals(data_host, indices_host, maxrowlength, dimension);

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
            double *vector_dev = NULL;
            double *data_dev = NULL;
            int *indices_dev = NULL;
            double *result = NULL;

            allocation(&vector_dev, dimension);
            allocation(&data_dev, dimension*maxrowlength);
            allocation(&indices_dev, dimension*maxrowlength);
            allocation(&result, dimension);

            copy_data(vector_host, vector_dev, dimension);
            copy_data(data_host, data_dev, dimension*maxrowlength);
            copy_data(indices_host, indices_dev, dimension*maxrowlength);

            invoke_gpu_overall(vector_dev, data_dev, indices_dev, result, dimension, maxrowlength);

            cleanup(vector_dev);
            cleanup(data_dev);
            cleanup(indices_dev);
            cleanup(result);
        }
        elapsed_overall = (timer_overall.stop()*1.0e3) / (float)iteration;
    }

//------------------------------------------------------------------------------------------------/
//                                   Zeitmessung Kernel
//------------------------------------------------------------------------------------------------/
    
    double *vector_dev = NULL;
    double *data_dev = NULL;
    int *indices_dev = NULL;
    double *result = NULL;

    allocation(&vector_dev, dimension);
    allocation(&data_dev, dimension*maxrowlength);
    allocation(&indices_dev, dimension*maxrowlength);
    allocation(&result, dimension);

    copy_data(vector_host, vector_dev, dimension);
    copy_data(data_host, data_dev, dimension*maxrowlength);
    copy_data(indices_host, indices_dev, dimension*maxrowlength);

    //=========================================//
    float elapsed_kernel = invoke_gpu_time(vector_dev, data_dev, indices_dev, result, dimension, maxrowlength, iteration);
    //>>>KERNEL<<<
    //=========================================//

    
    matvec_check_result(result, vector_host, data_host, indices_host, dimension, maxrowlength);
    cleanup(vector_dev);
    cleanup(data_dev);
    cleanup(indices_dev);
    cleanup(result);
 
//================================================================================================/
//                                         Evaluieren
//================================================================================================/

    double schalter = 0.0;
    performance(dimension, maxrowlength, elapsed_overall, elapsed_kernel, schalter, iteration);
  
    delete[] vector_host;
    delete[] data_host;
    delete[] indices_host;

    return 0;
}

