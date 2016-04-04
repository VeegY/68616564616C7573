#include "../src/include/fullvectorgpu.hpp"
#include "../src/include/fullvector.hpp"
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
int slicedgputest()
{
    srand (static_cast <unsigned> (time(0)));
    const size_t N=40;
    Icarus::FullVector<double> vec1(N), vec2(N), vec3(N);
    Icarus::FullVectorGpu<double> gvec1(N), gvec2(N), gvec3(N);
    size_t dimloc = vec1.get_dim();
    if (gvec1.get_dim()!=vec1.get_dim())
    {
        LOG_ERROR("get_dim failed");
    }
    double constdouble (static_cast <double> (rand()) / static_cast <double> (RAND_MAX));
    vec1.fill_const(constdouble);
    gvec1.fill_const(constdouble);
    for (size_t i(0); i<dimloc; i++)
    {
        double r = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
        vec2[i]= r;
        gvec2[i]= r;
    }
    Icarus::FullVectorGpu<double> gvec4(gvec1); //test copy constructor
    Icarus::FullVector<double> vec4(vec1); //test copy constructor
    gvec3=gvec1;
    vec3 =vec1;
    for (size_t i(0); i<dimloc; i++)
    {
        if (vec4[i]!=gvec4[i])
        {
             LOG_ERROR("copy-constructor failed");
        }
        if (vec3[i]!=gvec3[i])
        {
             LOG_ERROR("assignment operator failed");
        }
        if (vec1[i]!=gvec1[i])
        {
             LOG_ERROR("fill_const failed");
        }
    }

    //check artihmetic operations
    double randdouble= static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
    double checktol = std::numeric_limits<double>::epsilon()*100;
    vec2.scal(randdouble);
    gvec2.scal(randdouble);

    vec3.axpy(3.3, vec2);
    gvec3.axpy(3.3, gvec2);

    for (size_t i(0); i<dimloc; i++)
    {
        if ((vec2[i] - gvec2[i])>=(vec2[i] + gvec2[i])*checktol*10)
        {
             LOG_ERROR("scal");
        }
        if ((vec3[i] - gvec3[i])>=(vec3[i] + gvec3[i])*checktol*10)
        {
             LOG_ERROR("axpy");
        }

    }

LOG_INFO("arithmetics passed");
    vec4.copy(vec1);
    gvec4.copy(gvec1);
    for (size_t i(0); i<dimloc; i++)
    {
        if (vec4[i]!=gvec4[i]) LOG_ERROR ("copy failed");
    }
LOG_INFO("copy passed");

    if (vec1.maxnorm() - gvec1.maxnorm() >= checktol * (vec1.maxnorm() + gvec1.maxnorm()))
    {
        LOG_ERROR("maxnorm failed; cpu: ",vec1.maxnorm(), "  gpu: ", gvec1.maxnorm());
    }
LOG_INFO("maxnorm passed");
    if (vec1.l2norm2() - gvec1.l2norm2() >= checktol * (vec1.l2norm2() + gvec1.l2norm2()))
    {
        LOG_ERROR("L2norm2 failed; cpu: ",vec1.l2norm2(), "  gpu: ", gvec1.l2norm2());
    }
LOG_INFO("l2norm passed");
    if (vec1.l2norm() - gvec1.l2norm() >= checktol * (vec1.l2norm() + gvec1.l2norm()))
    {
        LOG_ERROR("L2norm failed; cpu: ",vec1.l2norm2(), "  gpu: ", gvec1.l2norm2());
    }

    return 0;
}

int main()
{
    return slicedgputest();
}
