#ifndef __MATHFUNCTION_HPP_
#define __MATHFUNCTION_HPP_

#include <limits>
#include <cassert>

namespace Icarus
{

class mathfunction
{
public:
    mathfunction(int type, double val=0.0):
        _type(type), _val(val) { }
    double eval(double x, double y, double z, int plane=0);
private:
    const int _type;
    const double _val;
};

}//namespace Icarus

#endif//__MATHFUNCTION_HPP_
