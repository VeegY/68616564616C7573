#ifndef __MATHFUNCTION_HPP_
#define __MATHFUNCTION_HPP_

namespace Icarus
{

class mathfunction
{
public:
    mathfunction(int type):
        _type(type) { }

    double eval(double x, double y, double z)
    {
        switch (_type)
        {
        case 0:
            return 0.0;
        case 1:
            return 1.0;
        case 2:
            return 3.14159265358979323846;
        default:
            return 0.0;
        }
        return 0.0;
    }
private:
    int _type;
};

}//namespace Icarus

#endif//__MATHFUNCTION_HPP_
