#include "include/assemblefem.hpp"

namespace Icarus
{

double assembleFem::getx(size_t index)
{
    int az = (index / Nx) / Ny;
    int ay = index / Nx - az*Ny;
    int ax = index - ay*Nx - az*Nx*Ny;
    return ((double) ax)*h;
}

double assembleFem::gety(size_t index)
{
    int az = (index / Nx) / Ny;
    int ay = index / Nx - az*Ny;
    return ((double)ay)*h;
}

double assembleFem::getz(size_t index)
{
    int az = (index / Nx) / Ny;
    return ((double) az)*h;
}

}//namespace Icarus
