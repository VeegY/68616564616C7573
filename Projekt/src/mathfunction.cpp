#include "include/mathfunction.hpp"

namespace Icarus
{

double mathfunction::eval(double x, double y, double z)
{
    switch (_type)
    {
    case 0: return _val;
    // 1: u, 2: f, 3: dirichlet, 4: neumann
    case 1: return (0.5-x)*(0.5-y)*(0.5-z);
    case 2: return 0.0;
    case 3: if (x < 1.0e-5 && x > -1e-5)
                return (0.5)*(0.5-y)*(0.5-z);
            if (x < 1.0+1e-5 && x > 1.0-1e-5)
                return (-0.5)*(0.5-y)*(0.5-z);
            if (y < 1.0e-5 && y > -1e-5)
                return (0.5-x)*(0.5)*(0.5-z);
            if (y < 1.0+1e-5 && y > 1.0-1e-5)
                return (0.5-x)*(-0.5)*(0.5-z);
            if (z < 1.0e-5 && z > -1e-5)
                return (0.5-x)*(0.5-y)*(0.5);
            if (z < 1.0+1e-5 && z > 1.0-1e-5)
                return (0.5-x)*(0.5-y)*(-0.5);
            assert(false || false); // nicht vorgesehen
    case 4: if (x < 1.0e-5 && x > -1e-5)
                return 0.5 - 0.5*(y+z) + y*z;
            if (x < 1.0+1e-5 && x > 1.0-1e-5)
                return - 0.5 + 0.5*(y+z) - y*z;
            if (y < 1.0e-5 && y > -1e-5)
                return 0.5 - 0.5*(x+z) + x*z;
            if (y < 1.0+1e-5 && y > 1.0-1e-5)
                return - 0.5 + 0.5*(x+z) - x*z;
            if (z < 1.0e-5 && z > -1e-5)
                return 0.5 - 0.5*(y+x) + y*x;
            if (z < 1.0+1e-5 && z > 1.0-1e-5)
                return - 0.5 + 0.5*(y+x) - y*x;
    // 5: u, 6: f, 7: dirichlet, 8: neumann
    case 5: return x*(1.0-x)*y*(1.0-y)*z*(1.0-z);
    case 6: return 2.0*(x*(1.0-x)*y*(1.0-y)+x*(1.0-x)*z*(1.0-z)+y*(1.0-y)*z*(1.0-z));
    case 7: return 0.0;
    case 8: assert(false); return 0.0; //TODO
    }
    assert(!true); // nicht vorgesehen
    return 0.0;
}

}//namespace Icarus
