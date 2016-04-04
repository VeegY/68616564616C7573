#include "../src/include/fullvectorgpu.hpp"

int main()
{
    Icarus::FullVectorGpu<double> testvec(10);
    LOG_INFO("passed");
    return 0;
}
