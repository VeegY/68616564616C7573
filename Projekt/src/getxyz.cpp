//TODO TODISCUSS ist es sinnvoll integer zu teilen um den ganzzahligen Quotiente zu bekommen?
//TODO TOCHECK zweite Variante wirklich nicht schneller?
#include "include/assemblefem.hpp"

namespace Icarus
{

double assembleFem::getx(size_t index)
{
    int az = (index / _nx) / _ny;
    int ay = index / _nx - az*_ny;
    int ax = index - ay*_nx - az*_nx*_ny;
    return ((double) ax)*_h;
//    return static_cast<double>(index - (index / _nx - (az)*_ny)*_nx - (az)*_nx*_ny)*_h;
}

double assembleFem::gety(size_t index)
{
    int az = (index / _nx) / _ny;
    int ay = index / _nx - az*_ny;
    return ((double)ay)*_h;
//    return static_cast<double>(index / _nx - ((index / _nx) / _ny)*_ny)*_h;
}

double assembleFem::getz(size_t index)
{
    int az = (index / _nx) / _ny;
    return ((double) az)*_h;
//    return static_cast<double>((index / _nx) / _ny)*_h;
}

}//namespace Icarus
