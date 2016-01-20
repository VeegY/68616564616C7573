#include <stdio.h>
#include <cmath>
#include <stdlib.h>
#include <iostream>

void random_ints(int *data,int *indices, int* fvec, int dim)
{
	for (int i = 0; i < dim;i++)
	{ 
  //Fuelle Reihe
        int entries_row=rand()%dim; //0-DIM Eintraege pro Reihe
        int place_current;
        int place_old=0;
        int range;
        for (int j = 0; j < dim; j++)
        {
	        //Fuelle i-te Stelle der Reihe
	        if (j < entries_row)
	        {
		        data[i*dim + j] = rand() % 1000;
		        range = dim - (entries_row - i) - place_old;
		        place_current = (rand() % range) + place_old + 1;
		        indices[i*dim + j] = place_current;
		        place_old = place_curent;
	        }
	        else
	        {
		        data[i*dim + j] = 0;
		        indices[i*dim + j] = 0;
	        }

        }
		fvec[i]=rand()%1000
    }
}

void set_values(int *datah, int *indicesh, int *fvech, int *datag, int *indicesg, int *fvecg, int dim)
{
	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			datag[i*dim + j] = datah[i*dim + j];
			indicesg[i*dim + j] = indicesh[i*dim + j];
		}
		fvecg[i] = fvech[i];
	}
}

void print_stuff(int *data, int *indices, int *fvec, int dim)
{
	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			printf("%i:%i - ", data[i*dim + j], indices[i*dim + j]);
		}
		printf(" --- fvec: %i\n", fvec[i]);
	}
}
