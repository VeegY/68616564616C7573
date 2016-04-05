#include "../src/include/slicedvectorgpu.hpp"
#include "../src/include/vtkwriter.hpp"
#include "../src/include/scalartraits.hpp"
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <limits>
/*
*test fuer slicedvectorgpu
*KEIN unit-test
*
*/
int slicedvectorgputest()
{
    srand (static_cast <unsigned> (time(0)));
    const size_t N=123481;
    Icarus::SlicedVectorGpu<double> vec1(N), vec2(N), vec4(N);
    Icarus::SlicedVectorGpu<double> vec6(N), vec7(N), vec8(N);
    size_t dimloc = vec1.get_dim_local();
    if (vec1.get_dim()!=N){
        LOG_ERROR("get_dim failed");
    }
    double constdouble (static_cast <double> (rand()) / static_cast <double> (RAND_MAX));
    vec4.fill_const(constdouble);
    for (size_t i(0); i<dimloc; i++)
    {
        double r = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
        int sign=2*(rand()%2)-1;
        vec1.set_local(i, r);
        vec6.set_local(i, sign*r);
    }
    Icarus::SlicedVectorGpu<double> vec3(vec1); //test copy constructor
    vec2=vec1;

    for (size_t i(0); i<dimloc; i++)
    {
        if (vec2.get_local(i)!=vec1.get_local(i))
        {
             LOG_ERROR("copy-constructor failed");
        }
        if (vec3.get_local(i)!=vec1.get_local(i))
        {
             LOG_ERROR("assignment operator failed");
        }
        if (vec4.get_local(i)!=constdouble)
        {
             LOG_ERROR("fill_const failed");
        }
    }

    //check artihmetic operations
    double randdouble= static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
    double maxnorm(0), l2norm2(0);
    double randdouble2= static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
    vec7=vec6;
    vec8=vec6;
    vec7.scal(randdouble);
    vec8.axpy(randdouble2, vec7);
    double checktol = std::numeric_limits<double>::epsilon();
    LOG_INFO("checktol: ", checktol);
    for (size_t i(0); i<dimloc; i++)
    {
        l2norm2+=vec6.get_local(i)*vec6.get_local(i);
        if (std::abs(vec6.get_local(i)) > maxnorm )
        {
            maxnorm=std::abs(vec6.get_local(i));
        }
        if (vec7.get_local(i)!=randdouble*vec6.get_local(i))
        {
            LOG_ERROR("scal failed ; value: ",vec7.get_local(i), "  reference value: ", randdouble*vec6.get_local(i));
        }
        if (vec8.get_local(i)!=vec6.get_local(i)+randdouble2*vec7.get_local(i))
        {
            LOG_INFO("axpy INFO; value: ",vec8.get_local(i), " difference: ", randdouble2*vec7.get_local(i)+vec6.get_local(i)-vec8.get_local(i));
        }
        if (std::abs(vec8.get_local(i) - (vec6.get_local(i)+randdouble2*vec7.get_local(i))) >=
             2*checktol) //*  std::abs(vec8.get_local(i) + (vec6.get_local(i)+randdouble2*vec7.get_local(i))))
        {
            LOG_ERROR("axpy failed; value: ",vec8.get_local(i), " difference: ", randdouble2*vec7.get_local(i)+vec6.get_local(i)-vec8.get_local(i));
        }
    }

    double maxnorm_glob, l2norm2_glob;
    MPI_SCALL(MPI_Allreduce(&maxnorm, &maxnorm_glob, 1, Icarus::ScalarTraits<double>::mpi_type, MPI_MAX, MPI_COMM_WORLD));
    MPI_SCALL(MPI_Allreduce(&l2norm2, &l2norm2_glob, 1, Icarus::ScalarTraits<double>::mpi_type, MPI_SUM, MPI_COMM_WORLD));
    if (maxnorm_glob!=vec6.maxnorm())
    {
        LOG_ERROR("maxnorm failed; value: ",vec6.maxnorm(), "  reference value: ", maxnorm_glob);
    }
    if (std::abs(l2norm2_glob-vec6.l2norm2())>=checktol*100*std::abs(l2norm2_glob+vec6.l2norm2()))
    {
        LOG_ERROR("L2norm2 failed; value: ",vec6.l2norm2(), "  reference value: ", l2norm2_glob, "  difference: ",l2norm2_glob-vec6.l2norm2());
    }
    double l2norm=std::sqrt(l2norm2_glob);
    if (std::abs(l2norm-vec6.l2norm())>=checktol*100*std::abs(l2norm+vec6.l2norm()))
    {
        LOG_ERROR("L2norm failed; value: ",vec6.l2norm(), "  reference value: ", l2norm, "  difference: ",l2norm-vec6.l2norm());
    }
    LOG_INFO("test passed");
    return 0;
}

int main()
{
    return slicedvectorgputest();
}
