#ifndef __MATHFUNCTION_HPP_
#define __MATHFUNCTION_HPP_

namespace Icarus
{

class mathfunction
{
public:
    mathfunction(int type, double val=0.0):
        _type(type), _val(val) { }

    double eval(double x, double y, double z)
    {
        switch (_type)
        {
        case 0: return _val;
        case 1: return 2.0*(x*(1.0-x)*y*(1.0-y)+x*(1.0-x)*z*(1.0-z)+y*(1.0-y)*z*(1.0-z));
        }
        return 0.0;
    }
private:
    const int _type;
    const double _val;
};

}//namespace Icarus

#endif//__MATHFUNCTION_HPP_
