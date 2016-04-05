#include "../src/include/fullvectorgpu.hpp"

int main()
{
    Icarus::SlicedVectorGpu<double> testvec(10000);
    LOG_INFO("passed");
    return 0;
}
