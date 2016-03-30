#include"include/timer.hpp"
#include <stdio.h>
#include <stdlib.h>
#include "mkl.h"

//------------------------------------------------------------------------------------------------/
//16384 - 65536 - 262144 - 1048576 - 4194304
#define dimension 262144
#define iterations 1000

void generate_stuff(int m, int k, double *A, int *idx, int *pntrb, int *pntre, double *x);

int main()
{
    Timer timer;
    int m(dimension), k(dimension);
    double *val = new double[7 * m];
    double *x = new double[k];
    double *y = new double[k];
    int *idx = new int[7 * m];
    int *pntrb = new int[m];
    int *pntre = new int[m];
    double alpha_v = 0.0;
    double beta_v = 0.0;
    double *alpha = &alpha_v;
    double *beta = &beta_v;
    

    generate_stuff(m, k, val, idx, pntrb, pntre, x);

    timer.start();
    mkl_dcsrmv('N', m, k, alpha, char *matdescra, val, idx, pntrb, pntre, x, beta, y);
    float elapsed = timer.stop()*1.0e3;

    return 0;
}


/*
void mkl_dcsrmv(char *transa, int *m, int *k, double *alpha, char *matdescra,
double *val, int *indx, int *pntrb, int *pntre, double *x, double *beta, double *y);

y := alpha*A*x + beta*y

transa = 'N' no transposed
m Number of rows matrix A
k Number of columns of the matrix A
matdescra = ??
val = array of non zero Elements of A(pref 0-7 per row for better comparing)
idx = column indices of those Elements; length = val;

pntrb, length m,
pntre, length m, pointer to 

*/

void generate_stuff(int m, int k, double *A,  int *idx, int *pntrb, int *pntre, double *x)
{
    srand(static_cast <unsigned> (time(0)));
    for (int i = 0; i < m; i++)
    {
        for (int n = 0; n < 7; n++)
        {
            A[m*7+n]= static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / 10));
            idx[m * 7 + n] = (rand() % (k - 1)) + 1;
        }
        pntrb[i] = 1 + 7 * i;
        pntre[i] = 8 + 7 * i;
    }
    for (int j = 0; i < k; i++)
    {
        x[j]= static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / 10));
    }
}

