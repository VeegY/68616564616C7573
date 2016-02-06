#include "../src/include/fullvector.hpp"
#include "../src/include/vtkwriter.hpp"
#include <cstdlib>
#include <ctime>
#include <cmath>
/*
*test fuer den vtk writer
*ausgabe im ordner out (muss existieren)
*KEIN unit-test
*
*/
int fullvectortest()
{
    srand (static_cast <unsigned> (time(0)));
    double r = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
	const size_t N=100000;
	const int randintmax;
	Icarus::FullVector<double> vec1(N), vec2(N), vec4(N);
	Icarus::FullVector<double> vecint1(N), vecint2(N), vecint3(N);
	if (vec1.get_dim())

	double constdouble (static_cast <double> (rand()) / static_cast <double> (RAND_MAX));
	vec4.fill_const(constdouble);
	for (size_t i(0); i<N; i++)
    {
        double r = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
        int sign=2*(rand()%2)-1;
        int randint=rand()%randintmax;
        vec1[i]=r;
        vecint1[i]=sign*randint;
    }
    FullVector<double> vec3(vec1); //test copy constructor
	vec2=vec1;
	for (size_t i(0); i<N; i++)
    {
        if (vec2[i]!=vec1[i])
        {
             LOG_ERROR("copy-constructor failed");
        }
        if (vec3[i]!=vec1[i])
        {
             LOG_ERROR("assignment operator failed");
        }
        if (vec4[i]!=constdouble)
        {
             LOG_ERROR("fill_const failed");
        }
    }

    //check artihmetic operations
    int randint=rand()%randintmax
    int maxnorm(0), l2norm2(0);
    vecint2=vecint1;
    vecint2.scal(randint);
    vecint3.axpy(randint, vec2);
    for (size_t i(0); i<N; i++)
    {
        l2norm2+=vecint1[i]*vecint1[i];
        if (abs(vecint1[i]) > maxnorm )
        {
            maxnorm=abs(vecint1[0]);
        }
        if (vecint2[i]!=randint*vecint1[i])
        {
             LOG_ERROR("scal failed");
        }
        if (vecint3[i]!=randint*vecint1[i]+vecint2[i])
        {
             LOG_ERROR("axpy failed");
        }
    }
    if (maxnorm!=vecint1.maxnorm())
    {
        LOG_ERROR("maxnorm failed");
    }
    if (l2norm2!=vecint1.l2norm2())
    {
        LOG_ERROR("L2norm2 failed");
    }
}

int main()
{
    int myrank;
    MPI_SCALL(MPI_Comm_rank(MPI_COMM_WORLD, &myrank));
    if (myrank == 0)
    {
        return fullvectortest();
    }
    return 0;
}
