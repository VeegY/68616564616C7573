namespace Icarus
{

class math_function
{
public:
    math_function(int type);
    double eval(double x, double y, double z);
private:
    int _type;
};

math_function::math_function(int type):
    _type(type)
{
}

double math_function::eval(double x, double y, double z)
{
    switch (number)
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

}//namespace Icarus
