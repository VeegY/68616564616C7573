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
                return -(0.5-x)*(0.5-y);
            if (plane == 2)
                return -(0.5-x)*(0.5-z);
            if (plane == 3)
                return -(0.5-y)*(0.5-z);
            assert(false); // nicht vorgesehen
    // 5: u, 6: f, 7: dirichlet, 8: neumann
    case 5: return x*(1.0-x)*y*(1.0-y)*z*(1.0-z);
    case 6: return 2.0*(x*(1.0-x)*y*(1.0-y)+x*(1.0-x)*z*(1.0-z)+y*(1.0-y)*z*(1.0-z));
    case 7: return 0.0;
    case 8: if (plane == 1)
                return x*(1.0-x)*y*(1.0-y)*(1.0-2.0*z);
            if (plane == 2)
                return x*(1.0-x)*(1.0-2.0*y)*z*(1.0-z);
            if (plane == 3)
                return (1.0-2.0*x)*y*(1.0-y)*z*(1.0-z);
            assert(false); // nicht vorgesehen
    // 9: u, 10: f, 11: dirichlet, 12: neumann
    case 9: return x;
    case 10:return 0.0;
    case 11:return x;
    case 12:if (plane == 1)
                return 0.0;
            if (plane == 2)
                return 0.0;
            if (plane == 3)
                return 1.0;
            assert(false); // nicht vorgesehen
    // 13: u, 14: f, 15: dirichlet, 16: neumann - femflow_test
    case 13:return 0.0;
    case 14:return 0.0;
    case 15:return 0.0;
    case 16:if (plane == 1)
                return 0.0;
            if (plane == 2)
                return 0.0;
            if (plane == 3)
                return 10.0;
            assert(false); // nicht vorgesehen
    case 99:if (plane == 1)
            {
                // Mitte der Luefters:
                // x = (0.4 + 0.633 ) / 2
                // y = (0.30147 + 0.53733) / 2
                // z = 0.12623
                // Aufbau: Deckel | prop | Deckel | Mitte | Deckel | prop | Deckel
                if ((x-0.5165)*(x-0.5165) + (y-0.4194)*(y-0.4194) < 0.117   // nah genug dran
                    && (x-0.5165)*(x-0.5165) + (y-0.4194)*(y-0.4194) > 0.5  // weit genug weg
                    && z > 0.12623-1.0e-3 && z < 0.12623+1.0e-3)            // Hoehe richtig
                    return 50.0;
                return 0.0; // nicht Luefter
            }
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
