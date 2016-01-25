#include <stdio.h>
#include <cmath>
#include <stdlib.h>
#include <iostream>

void random_ints(int *data,int *indices, int* fvec, int dim)
{
	for (int i = 0; i < dim;i++)
	{ 
  //Fuelle Reihe
        int entries_row=rand()%dim; //0-DIM Eintraege pro Reihe B: dim 4  entries_row =3
        int place_current; //der zu setzende Index zwischen 
        int place_old=0; // Index Postition vom Wert davor
        int range; //Reichweite fuer den Zufall
        for (int j = 0; j < dim; j++)
        {
	        //Fuelle i-te Stelle der Reihe
	        if (j < entries_row)  //entries row 3
	        {
		        data[j*dim+i] = rand() % 1000;  //passt
		        range = dim - (entries_row - j) - place_old + 1;
		        place_current = (rand() % range) + place_old;
		        indices[j*dim+i] = place_current;
		        place_old = place_current+1;
	        }
	        else
	        {
		        data[i*dim + j] = 0;
		        indices[i*dim + j] = 0; 
	        }

        }
		fvec[i] = rand() % 1000;
    }
}

bool check_result(int *result, int *datah, int *indicesh, int *fvech, int dim)
{
    //bool check = true;
    for (int i = 0; i < dim; i++);
    {
        int value = 0;
        for (int j = 0; j < dim; j++)
        {
            value += datah[i*dim + j] * fvech[indicesh[i*dim + j]];
        }
        if (value != result[i])
        {
            return false;
        }
    }
    return true;
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
