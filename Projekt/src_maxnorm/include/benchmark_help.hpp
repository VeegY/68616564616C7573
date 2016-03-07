#include <stdio.h>
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <ctime>
#include <cstdlib>
#define WARNING "\e[30;41;1m"
#define GREEN "\e[32;1m"
#define RESET "\e[0m"

void vec_float_one(float *vecone, int dim_local)
{
    srand(static_cast <unsigned> (time(0)));
    float value_one;
    for (int k = 0; k < dim_local; k++)
    {
        value_one = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / 10));
        vecone[k] = value_one;
    }

}



void check_result_maxnorm(float *result, float *vector, int dim_local, char a)
{
    float diff, value = 0.0;
    bool check = true;
    for (int k = 0; k < dim_local; k++)
    {
        compare = vector[k];
        if (compare < 0)
        {
            compare = -compare;
        }
        
        if (compare > value)
        {
            value = compare;
        }
    }
    diff = value - result[0];
    printf(GREEN "DIFF: %f\nVALUE: %f\nRESULT: %f\n" RESET, diff, value, result[0]);
    if (diff > 1.0e-2 || diff < -1.0e-2)
    {
        check = false;
    }

    if (check)
    {
        printf(GREEN "%c_Kernel outcome true\n" RESET, a);
    }
    else printf(WARNING "%c_Kernel outcome false\n", a);
}



void set_values_maxnorm(float *veconeh, float *veconeg, int dim_local)
{
    for (int k = 0; k < dim_local; k++)
    {
        veconeg[k] = veconeh[k];
    }
}

void print_vec(float *vecone, float *vectwo, int dim_local)
{
    for (int k = 0; k < dim_local; k++)
    {
        printf("Vector One: %f ~~ Vector Two: %f\n", vecone[k], vectwo[k]);
    }
}

void print_vec_one(float *vecone, int dim_local)
{
    for (int k = 0; k < dim_local; k++)
    {
        printf("Vector One: %f\n", vecone[k]);
    }
}


void print_time(float *ukt, float *uot, float *zkt, float *zot,int runs)
{
    float uktime=0, uotime=0, zktime=0, zotime=0;
    for(int n=0;n<runs;n++)
    {
      uktime += ukt[n];
      uotime += uot[n];
      zktime += zkt[n];
      zotime += zot[n];
    }
    uktime = (uktime/runs)*1000;
    uotime = (uotime/runs)*1000;
    zktime = (zktime/runs)*1000;
    zotime = (zotime/runs)*1000;

    printf("UK: %fms - UO: %fms - ZK: %fms - ZO: %fms\n",uktime,uotime,zktime,zotime);
}

