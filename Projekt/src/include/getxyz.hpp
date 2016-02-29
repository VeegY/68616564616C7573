#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <cmath>
namespace Icarus
{

<<<<<<< HEAD

    double getx(size_t index,
		double h, size_t nx, size_t ny)
=======
int getx(size_t index, size_t nx, size_t ny)
>>>>>>> 56616b8691e7ea15852f116f5c07e7949baa1f8b
	{
	int	az = (index / nx) / ny;
	int	ay = index / nx - az*ny;
	int	ax = index - ay*nx - az*nx*ny;

<<<<<<< HEAD
    return ((double) ax)*h;
	}

	double gety(size_t index,
		double h, size_t nx, size_t ny)
=======
    	return ax;
	}

int gety(size_t index, size_t nx, size_t ny)
>>>>>>> 56616b8691e7ea15852f116f5c07e7949baa1f8b
	{
	int	az = (index / nx) / ny;
	int	ay = index / nx - az*ny;

<<<<<<< HEAD
    return ((double)ay)*h;
	}

	double getz(size_t index,
		double h, size_t nx, size_t ny)
=======
    	return ay;
	}

int getz(size_t index, size_t nx, size_t ny)
>>>>>>> 56616b8691e7ea15852f116f5c07e7949baa1f8b
	{
    	int az = (index / nx) / ny;

<<<<<<< HEAD
    return ((double) az)*h;
=======
    	return az;
>>>>>>> 56616b8691e7ea15852f116f5c07e7949baa1f8b
	}
}
