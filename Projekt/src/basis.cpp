//TODO TODISCUSS:
// immer direkt return ...; ? also zwsp und x0,y0,z0 umgehen

#include "include/assemblefem.hpp"

namespace Icarus
{

// Basis auf Referenzelement [0,h]^3 evaluiert in den Gausspunkten
std::vector<double> assembleFem::evaluated_Basis3d(int A)
{
    assert(A >= 0 && A < 8);

    // Fuer Quadraturpunkte
    std::vector<double> zwsp(27);
//    std::vector<double> X(27), Y(27), Z(27), zwsp(27);
//    X = get_quadrature_xpoints();
//    Y = get_quadrature_ypoints();
//    Z = get_quadrature_zpoints();

    switch(A)
    {
    case 0: for(int i(0); i < 27; ++i)
                zwsp[i]= ((_h - _quadpoints_3d_x[i]) * (_h - _quadpoints_3d_y[i]) * (_h - _quadpoints_3d_z[i])) / (_h*_h*_h);
            break;
    case 1: for(int i(0); i < 27; ++i)
                zwsp[i]= (_quadpoints_3d_x[i] * (_h - _quadpoints_3d_y[i]) * (_h - _quadpoints_3d_z[i])) / (_h*_h*_h);
            break;
    case 2: for(int i(0); i < 27; ++i)
                zwsp[i]= ((_h - _quadpoints_3d_x[i])  * _quadpoints_3d_y[i]  * (_h - _quadpoints_3d_z[i])) / (_h*_h*_h);
            break;
    case 3: for(int i(0); i < 27; ++i)
                zwsp[i]= (_quadpoints_3d_x[i] * _quadpoints_3d_y[i] * (_h - _quadpoints_3d_z[i])) / (_h*_h*_h);
            break;
    case 4: for(int i(0); i < 27; ++i)
                zwsp[i]= ((_h - _quadpoints_3d_x[i]) * (_h - _quadpoints_3d_y[i]) * _quadpoints_3d_z[i]) / (_h*_h*_h);
            break;
    case 5: for(int i(0); i < 27; ++i)
                zwsp[i]= (_quadpoints_3d_x[i] * (_h - _quadpoints_3d_y[i]) * _quadpoints_3d_z[i]) / (_h*_h*_h);
            break;
    case 6: for(int i(0); i < 27; ++i)
                zwsp[i]= ((_h - _quadpoints_3d_x[i]) * _quadpoints_3d_y[i] * _quadpoints_3d_z[i]) / (_h*_h*_h);
            break;
    case 7: for(int i(0); i < 27; ++i)
                zwsp[i]= (_quadpoints_3d_x[i] * _quadpoints_3d_y[i] * _quadpoints_3d_z[i]) / (_h*_h*_h);
            break;
    }

    return zwsp;
}

// Gradienten der Basis auf Referenzelement [0,h]^3 evaluiert in den Gausspunkten
// ALTERNATIV: Matrix mit Dimension 3x27
std::vector<std::vector<double>> assembleFem::evaluated_gradient_Basis3d(int A)
{
    assert(A >= 0 && A < 8);

    // Fuer Quadraturpunkte
    std::vector<double> zwsp1(27), zwsp2(27), zwsp3(27);
//    std::vector<double> X(27), Y(27), Z(27), zwsp1(27), zwsp2(27), zwsp3(27);
//    X = get_quadrature_xpoints();
//    Y = get_quadrature_ypoints();
//    Z = get_quadrature_zpoints();

    //TODO _h*_h*_h ersetzen
    switch(A)
    {
    case 0: for(int i(0); i < 27; ++i)
            {
                zwsp1[i] = - ((_h - _quadpoints_3d_y[i]) * (_h - _quadpoints_3d_z[i])) / (_h*_h*_h);
                zwsp2[i] = - ((_h - _quadpoints_3d_x[i]) * (_h - _quadpoints_3d_z[i])) / (_h*_h*_h);
                zwsp3[i] = - ((_h - _quadpoints_3d_x[i]) * (_h - _quadpoints_3d_y[i])) / (_h*_h*_h);
            }
            break;
    case 1: for(int i(0); i < 27; ++i)
            {
                zwsp1[i] = (_h - _quadpoints_3d_y[i]) * (_h - _quadpoints_3d_z[i]) / (_h*_h*_h);
                zwsp2[i] = - (_quadpoints_3d_x[i] * (_h - _quadpoints_3d_z[i])) / (_h*_h*_h);
                zwsp3[i] = - (_quadpoints_3d_x[i] * (_h - _quadpoints_3d_y[i])) / (_h*_h*_h);
            }
            break;
    case 2: for(int i(0); i < 27; ++i)
            {
                zwsp1[i] = - (_quadpoints_3d_y[i] * (_h - _quadpoints_3d_z[i])) / (_h*_h*_h);
                zwsp2[i] = ((_h - _quadpoints_3d_x[i]) * (_h - _quadpoints_3d_z[i])) / (_h*_h*_h);
                zwsp3[i] = - ((_h - _quadpoints_3d_x[i]) * _quadpoints_3d_y[i]) / (_h*_h*_h);
            }
            break;
    case 3: for(int i(0); i < 27; ++i)
            {
                zwsp1[i] = (_quadpoints_3d_y[i] * (_h - _quadpoints_3d_z[i])) / (_h*_h*_h);
                zwsp2[i] = (_quadpoints_3d_x[i] * (_h - _quadpoints_3d_z[i])) / (_h*_h*_h);
                zwsp3[i] = - (_quadpoints_3d_x[i] * _quadpoints_3d_y[i]) / (_h*_h*_h);
            }
            break;
    case 4: for(int i(0); i < 27; ++i)
            {
                zwsp1[i] = - ((_h - _quadpoints_3d_y[i]) * _quadpoints_3d_z[i]) / (_h*_h*_h);
                zwsp2[i] = - ((_h - _quadpoints_3d_x[i]) * _quadpoints_3d_z[i]) / (_h*_h*_h);
                zwsp3[i] = ((_h - _quadpoints_3d_x[i]) * (_h - _quadpoints_3d_y[i])) / (_h*_h*_h);
            }
            break;
    case 5: for(int i(0); i < 27; ++i)
            {
                zwsp1[i] = ((_h - _quadpoints_3d_y[i]) * _quadpoints_3d_z[i]) / (_h*_h*_h);
                zwsp2[i] = - (_quadpoints_3d_x[i] * _quadpoints_3d_z[i]) / (_h*_h*_h);
                zwsp3[i] = (_quadpoints_3d_x[i] * (_h - _quadpoints_3d_y[i])) / (_h*_h*_h);
            }
            break;
    case 6: for(int i(0); i < 27; ++i)
            {
                zwsp1[i] = - (_quadpoints_3d_y[i] * _quadpoints_3d_z[i]) / (_h*_h*_h);
                zwsp2[i] = ((_h - _quadpoints_3d_x[i]) * _quadpoints_3d_z[i]) / (_h*_h*_h);
                zwsp3[i] = ((_h - _quadpoints_3d_x[i]) * _quadpoints_3d_y[i]) / (_h*_h*_h);
            }
            break;
    case 7: for(int i(0); i < 27; ++i)
            {
                zwsp1[i] = (_quadpoints_3d_y[i] * _quadpoints_3d_z[i]) / (_h*_h*_h);
                zwsp2[i] = (_quadpoints_3d_x[i] * _quadpoints_3d_z[i]) / (_h*_h*_h);
                zwsp3[i] = (_quadpoints_3d_x[i] * _quadpoints_3d_y[i]) / (_h*_h*_h);
            }
            break;
    }

    return {zwsp1, zwsp2, zwsp3};
}


// Basis auf Referenzelement [0,h]^2 evaluiert in den Gausspunkten
std::vector<double> assembleFem::evaluated_Basis2d_1(int A)
{
    assert(A >= 0 && A < 4);

    // Fuer Quadraturpunkte
    std::vector<double> zwsp(9);
//    std::vector<double> R1(9), R2(9);
//    R1 = get_quadrature_xpoints_2d_1();
//    R2 = get_quadrature_ypoints_2d_1();

    switch (A)
    {
    case 0: for(int i(0); i < 9; ++i)
                zwsp[i] = ((_h - _quadpoints_2d_1[i]) * (_h - _quadpoints_2d_2[i])) / (_h*_h);
            break;
    case 1: for(int i(0); i < 9; ++i)
                zwsp[i] = (_quadpoints_2d_1[i] * (_h - _quadpoints_2d_2[i])) / (_h*_h);
            break;
    case 2: for(int i(0); i < 9; ++i)
                zwsp[i] = ((_h - _quadpoints_2d_1[i]) * _quadpoints_2d_2[i]) / (_h*_h);
            break;
    case 3: for(int i(0); i < 9; ++i)
                zwsp[i] = (_quadpoints_2d_1[i] * _quadpoints_2d_2[i]) / (_h*_h);
            break;
    }

    return zwsp;
}

// Basis auf Referenzelement [0,h]^2 evaluiert in den Gausspunkten
std::vector<double> assembleFem::evaluated_Basis2d_2(int A)
{
    assert(A == 0 || A == 1 || A == 4 || A == 5);

    // Fuer Quadraturpunkte
    std::vector<double> zwsp(9);
//    std::vector<double> R1(9), R2(9);
//    R1 = get_quadrature_xpoints_2d_2();
//    R2 = get_quadrature_zpoints_2d_2();

    switch (A)
    {
    case 0: for(int i(0); i < 9; ++i)
                zwsp[i] = ((_h - _quadpoints_2d_1[i]) * (_h - _quadpoints_2d_2[i])) / (_h*_h);
            break;
    case 1: for(int i(0); i < 9; ++i)
                zwsp[i] = (_quadpoints_2d_1[i] * (_h - _quadpoints_2d_2[i])) / (_h*_h);
            break;
    case 4: for(int i(0); i < 9; ++i)
                zwsp[i] = ((_h - _quadpoints_2d_1[i]) * _quadpoints_2d_2[i]) / (_h*_h);
            break;
    case 5: for(int i(0); i < 9; ++i)
               zwsp[i] = (_quadpoints_2d_1[i] * _quadpoints_2d_2[i]) / (_h*_h);
            break;
    }

    return zwsp;
}

std::vector<double> assembleFem::evaluated_Basis2d_3(int A)
{
    assert(A == 0 || A == 2 || A == 4 || A == 6);

    // Fuer Quadraturpunkte
    std::vector<double> zwsp(9);
//    std::vector<double> R1(9), R2(9);
//    R1 = get_quadrature_ypoints_2d_3();
//    R2 = get_quadrature_zpoints_2d_3();

    switch (A)
    {
    case 0: for(int i(0); i < 9; ++i)
                zwsp[i] = ((_h - _quadpoints_2d_1[i]) * (_h - _quadpoints_2d_2[i])) / (_h*_h);
            break;
    case 2: for(int i(0); i < 9; ++i)
                zwsp[i] = (_quadpoints_2d_1[i] * (_h - _quadpoints_2d_2[i])) / (_h*_h);
            break;
    case 4: for(int i(0); i < 9; ++i)
                zwsp[i] = ((_h - _quadpoints_2d_1[i]) * _quadpoints_2d_2[i]) / (_h*_h);
            break;
    case 6: for(int i(0); i < 9; ++i)
               zwsp[i] = (_quadpoints_2d_1[i] * _quadpoints_2d_2[i]) / (_h*_h);
            break;
    }

    return zwsp;
}

}//namespace Icarus
