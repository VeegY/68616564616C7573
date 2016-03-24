#include "include/mathfunction.hpp"

namespace Icarus
{

double mathfunction::eval(double x, double y, double z, int plane)
{
    switch (_type)
    {
    case 0: return _val;
    // 1: u, 2: f, 3: dirichlet, 4: neumann
    case 1: return (0.5-x)*(0.5-y)*(0.5-z);
    case 2: return 0.0;
    case 3: return (0.5-x)*(0.5-y)*(0.5-z);
            if (x < 1.0e-5 && x > -1.0e-5)
                return (0.5)*(0.5-y)*(0.5-z);
            if (x < 1.0+1.0e-5 && x > 1.0-1.0e-5)
                return (-0.5)*(0.5-y)*(0.5-z);
            if (y < 1.0e-5 && y > -1.0e-5)
                return (0.5-x)*(0.5)*(0.5-z);
            if (y < 1.0+1.0e-5 && y > 1.0-1.0e-5)
                return (0.5-x)*(-0.5)*(0.5-z);
            if (z < 1.0e-5 && z > -1.0e-5)
                return (0.5-x)*(0.5-y)*(0.5);
            if (z < 1.0+1.0e-5 && z > 1.0-1.0e-5)
                return (0.5-x)*(0.5-y)*(-0.5);
            assert(false); // nicht vorgesehen
    case 4: if (plane == 1)
                return -(0.5-y)*(0.5-z);
            if (plane == 2)
                return -(0.5-x)*(0.5-z);
            if (plane == 3)
                return -(0.5-x)*(0.5-y);
            assert(false); // nicht vorgesehen
    // 5: u, 6: f, 7: dirichlet, 8: neumann
    case 5: return x*(1.0-x)*y*(1.0-y)*z*(1.0-z);
    case 6: return -2.0*(x*(1.0-x)*y*(1.0-y)+x*(1.0-x)*z*(1.0-z)+y*(1.0-y)*z*(1.0-z));
    case 7: return 0.0;
    case 8: if (plane == 1)
                return (1.0-2.0*x)*y*(1.0-y)*z*(1.0-z);
            if (plane == 2)
                return x*(1.0-x)*(1.0-2.0*y)*z*(1.0-z);
            if (plane == 3)
                return x*(1.0-x)*y*(1.0-y)*(1.0-2.0*z);
            assert(false); // nicht vorgesehen
    // 9: u, 10: f, 11: dirichlet, 12: neumann
    case 9: return x*x;
    case 10:return 2.0;
    case 11:return x*x;
    case 12:if (plane == 1)
                return 2*x;
            if (plane == 2)
                return 0.0;
            if (plane == 3)
                return 0.0;
            assert(false); // nicht vorgesehen
    }
    assert(!true); // nicht vorgesehen
    return 0.0;
}

}//namespace Icarus
