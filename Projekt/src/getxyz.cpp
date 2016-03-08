#include "include/assemblefem.hpp"

namespace Icarus
{

double assembleFem::getx(size_t index)
{
    int az = (index / _nx) / _ny;
    int ay = index / _nx - az*_ny;
    int ax = index - ay*_nx - az*_nx*_ny;
    return ((double) ax)*_h;
}

double assembleFem::gety(size_t index)
{
    int az = (index / _nx) / _ny;
    int ay = index / _nx - az*_ny;
    return ((double)ay)*_h;
}

double assembleFem::getz(size_t index)
{
    int az = (index / _nx) / _ny;
    return ((double) az)*_h;
}

}//namespace Icarus
