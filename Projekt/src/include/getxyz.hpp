#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <cmath>
namespace Icarus
{


    double getx(size_t index,
		double h, size_t nx, size_t ny)
	{
	int	az = (index / nx) / ny;
	int	ay = index / nx - az*ny;
	int	ax = index - ay*nx - az*nx*ny;

    return ((double) ax)*h;
	}

	double gety(size_t index,
		double h, size_t nx, size_t ny)
	{
	int	az = (index / nx) / ny;
	int	ay = index / nx - az*ny;

    return ((double)ay)*h;
	}

	double getz(size_t index,
		double h, size_t nx, size_t ny)
	{
    int az = (index / nx) / ny;

    return ((double) az)*h;
	}
}
