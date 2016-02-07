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
	const size_t N=100000;
	Icarus::FullVector<double> vec1(N), vec2(N), vec4(N);
	Icarus::FullVector<double> vec6(N), vec7(N), vec8(N);
	if (vec1.get_dim()!=N){
        LOG_ERROR("get_dim failed");
	}
	double constdouble (static_cast <double> (rand()) / static_cast <double> (RAND_MAX));
	vec4.fill_const(constdouble);
	for (size_t i(0); i<N; i++)
    {
        double r = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
        int sign=2*(rand()%2)-1;
        vec1[i]=r;
        vec6[i]=sign*r;
    }
    Icarus::FullVector<double> vec3(vec1); //test copy constructor
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
    double randdouble= static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
    double maxnorm(0), l2norm2(0);
    vec7=vec6;
    vec8=vec6;
    vec7.scal(randdouble);
    vec8.axpy(randdouble, vec7);
    for (size_t i(0); i<N; i++)
    {
        l2norm2+=vec6[i]*vec6[i];
        if (abs(vec6[i]) > maxnorm )
        {
            maxnorm=std::abs(vec6[i]);
        }
        if (vec7[i]!=randdouble*vec6[i])
        {
             LOG_DEBUG("scal failed ; value: ",vec7[i], "  reference value: ", randdouble*vec6[i]);
        }
        if (vec8[i]!=vec6[i]+randdouble*vec7[i])
        {
             LOG_DEBUG("axpy failed; value: ",vec8[i], "  reference value: ", randdouble*vec6[i]+vec7[i]);
        }
    }
    if (maxnorm!=vec6.maxnorm())
    {
        LOG_DEBUG("maxnorm failed; value: ",vec6.maxnorm(), "  refernce value: ", maxnorm);
    }
    if (l2norm2!=vec6.l2norm2())
    {
        LOG_DEBUG("L2norm2 failed; value: ",vec6.l2norm2(), "  refernce value: ", l2norm2);
    }
    return 0;
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
