#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <cmath>
namespace Icarus;
{


    int getx(size_t index,
		size_t nx, size_t ny)
	{
	int	az = (index / nx) / ny;
	int	ay = index / nx - az*ny;
	int	ax = index - ay*nx - az*nx*ny;

    return ax;
	}

	int gety(size_t index,
		size_t nx, size_t ny)
	{
	int	az = (index / nx) / ny;
	int	ay = index / nx - az*ny;

    return ay;
	}

	int getz(size_t index,
		size_t nx, size_t ny)
	{
    int az = (index / nx) / ny;

    return az;
	}
}
