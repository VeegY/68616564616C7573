#include <stdio.h>
#include <cmath>
#include <stdlib.h>
void random_ints(int *data,int *indices, int dim)
{
  //Fuelle Reihe
  int entries_row=rand()%dim; //0-DIM Eintraege pro Reihe
  int place_current;
  int place_old=0;
  int range;
  for(int i=0;i<dim;i++)
  {
    //Fuelle i-te Stelle der Reihe
    if(i<entries_row)
    {
      data[i]=rand();
      range=dim-(entries_row-i)-place_old;
      place_current=(rand()%range)+place_old+1;
      indices[i]=place_current;
      place_old=place_curent;
    }
    else
    {
      data[i]=0;
      indices[i]=0;
    }

    //Indices
  }
}

void print_stuff(int *data, int *indices, int *fvec, int dim)
{
}
