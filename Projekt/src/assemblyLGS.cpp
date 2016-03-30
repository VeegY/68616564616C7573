#include "include/assemblefem.hpp"
#include "include/mathfunction.hpp"

#include "include/logger.hpp"
#include <iostream>

namespace Icarus
{

void assembleFem::assemble(DistEllpackMatrix<double>& Matrix, SlicedVector<double>& rhs,
    std::vector<char>& disc_points, mathfunction f, mathfunction g, mathfunction h)
{
    //TODO alles in size_t umwandeln
    int fron(Matrix.first_row_on_node()), lron(fron + Matrix.get_dim_local()-1);

    LOG_INFO(disc_points.size(), ", ", _nx*_ny*_nz);
    //TODO: vorlaeufig, wieder loeschen
    bool Dirichlet(false);
    bool Neumann(true);
    //TODO: vorlaeufig, wieder loeschen

    Matrix.prepare_sequential_fill(27);
    int rowlength(0);

    int Zeile;
    std::vector<double> RHS(_nx*_ny*_nz);

LOG_INFO("assembled 0%");
    //Ecke 1
    _e.clear(); _e.resize(1);
    _A.clear(); _A.resize(1);
    rowlength = 8;
    Zeile=0;
if (Zeile >= fron && Zeile <= lron)
{
//    if(Dirichlet)
    {
        Matrix.sequential_fill(Zeile, 1.0);
        Matrix.end_of_row();
        RHS[Zeile]= g.eval(getx(Zeile), gety(Zeile), getz(Zeile));
    }
//    else
//    {
//        _e[0]=Zeile; _A[0]=0;
//        assemblyMatrixRow(rowlength);
//        for (int m(0); m < rowlength; ++m)
//            Matrix.sequential_fill(_column[m], _value[m]);
//        Matrix.end_of_row();
//        RHS[Zeile] = assemblyRHSLoad(f);
//        if(Neumann)
//        {
//            RHS[Zeile] += assemblyRHSNeumann(1, false, h);
//            RHS[Zeile] += assemblyRHSNeumann(2, false, h);
//            RHS[Zeile] += assemblyRHSNeumann(3, false, h);
//        }
//    }
}

    //Kante 1:
    _e.resize(2);
    _A.resize(2);
    rowlength = 12;
    for(int i(1); i<_nx-1;i++)
    {
        Zeile++;
if (Zeile >= fron && Zeile <= lron)
{
        if(Dirichlet)
        {
            Matrix.sequential_fill(Zeile, 1.0);
            Matrix.end_of_row();
            RHS[Zeile] = g.eval(getx(Zeile), gety(Zeile), getz(Zeile));
        }
        else
        {
            _e[0]=Zeile-1; _A[0]=1;
            _e[1]=Zeile; _A[1]=0;
            assemblyMatrixRow(rowlength);
            for (int m(0); m < rowlength; ++m)
                Matrix.sequential_fill(_column[m], _value[m]);
            Matrix.end_of_row();
            RHS[Zeile] = assemblyRHSLoad(f);
            if(Neumann)
            {
                RHS[Zeile] += assemblyRHSNeumann(1, false, h);
                RHS[Zeile] += assemblyRHSNeumann(2, false, h);
            }
        }
}
    }//close I-Schleife (X-Achse)

    //Ecke 2
    _e.resize(1);
    _A.resize(1);
    rowlength = 8;
    Zeile++;
if (Zeile >= fron && Zeile <= lron)
{ //Zeile sollte hier y-1 sein
    if(Dirichlet)
    {
        Matrix.sequential_fill(Zeile, 1.0);
        Matrix.end_of_row();
        RHS[Zeile]= g.eval(getx(Zeile), gety(Zeile), getz(Zeile));
    }
    else
    {
        _e[0]=Zeile-1; _A[0]=1;
        assemblyMatrixRow(rowlength);
        for (int m(0); m < rowlength; ++m)
            Matrix.sequential_fill(_column[m], _value[m]);
        Matrix.end_of_row();
        RHS[Zeile] = assemblyRHSLoad(f);
        if(Neumann)
        {
            RHS[Zeile] += assemblyRHSNeumann(1, false, h);
            RHS[Zeile] += assemblyRHSNeumann(2, false, h);
            RHS[Zeile] += assemblyRHSNeumann(3, true, h);
        }
    }
}
    for(int j(1); j<_ny-1;j++)
    {
        //Kante 5
        _e.resize(2);
        _A.resize(2);
        rowlength = 12;
        Zeile++;
if (Zeile >= fron && Zeile <= lron)
{
        if(Dirichlet)
        {
            Matrix.sequential_fill(Zeile, 1.0);
            Matrix.end_of_row();
            RHS[Zeile]= g.eval(getx(Zeile), gety(Zeile), getz(Zeile));
        }
        else
        {
            _e[0]=Zeile-y; _A[0]=2;
            _e[1]=Zeile; _A[1]=0;
            assemblyMatrixRow(rowlength);
            for (int m(0); m < rowlength; ++m)
                Matrix.sequential_fill(_column[m], _value[m]);
            Matrix.end_of_row();
            RHS[Zeile] = assemblyRHSLoad(f);
            if(Neumann)
            {
                RHS[Zeile] += assemblyRHSNeumann(1, false, h);
                RHS[Zeile] += assemblyRHSNeumann(3, false, h);
            }
        }
}
        //Flaeche 1
        _e.resize(4);
        _A.resize(4);
        rowlength = 18;
        for(int i(1); i<_nx-1;i++)
        {
            Zeile++;
if (Zeile >= fron && Zeile <= lron)
{
            if(Dirichlet)
            {
                Matrix.sequential_fill(Zeile, 1.0);
                Matrix.end_of_row();
                RHS[Zeile]= g.eval(getx(Zeile), gety(Zeile), getz(Zeile));
            }
            else
            {
                _e[0]=Zeile -y-1; _A[0]=3;
                _e[1]=Zeile -y; _A[1]=2;
                _e[2]=Zeile -1; _A[2]=1;
                _e[3]=Zeile; _A[3]=0;
                assemblyMatrixRow(rowlength);
                for (int m(0); m < rowlength; ++m)
                    Matrix.sequential_fill(_column[m], _value[m]);
                Matrix.end_of_row();
                RHS[Zeile] = assemblyRHSLoad(f);
                if(Neumann)
                {
                    RHS[Zeile] += assemblyRHSNeumann(1, false, h);
                }
            }
}
        } //close I-Schleife (X-Achse)

        //Kante: 6
        _e.resize(2);
        _A.resize(2);
        rowlength = 12;
        Zeile++;
if (Zeile >= fron && Zeile <= lron)
{
        if(Dirichlet)
        {
            Matrix.sequential_fill(Zeile, 1.0);
            Matrix.end_of_row();
            RHS[Zeile]= g.eval(getx(Zeile), gety(Zeile), getz(Zeile));
        }
        else
        {
            _e[0]=Zeile -1 -y; _A[0]=3;
            _e[1]=Zeile -1; _A[1]=1;
            assemblyMatrixRow(rowlength);
            for (int m(0); m < rowlength; ++m)
                Matrix.sequential_fill(_column[m], _value[m]);
            Matrix.end_of_row();
            RHS[Zeile] = assemblyRHSLoad(f);
            if(Neumann)
            {
                RHS[Zeile] += assemblyRHSNeumann(1, false, h);
                RHS[Zeile] += assemblyRHSNeumann(3, true, h);
            }
        }
}
    } //close J-Schleife (Y-Achse)

    //Ecke 3:
    _e.resize(1);
    _A.resize(1);
    rowlength = 8;
    Zeile++;
if (Zeile >= fron && Zeile <= lron)
{ //Zeile sollte hier (_ny-1)*y sein
    if(Dirichlet)
    {
        Matrix.sequential_fill(Zeile, 1.0);
        Matrix.end_of_row();
        RHS[Zeile]= g.eval(getx(Zeile), gety(Zeile), getz(Zeile));
    }
    else
    {
        _e[0]=Zeile-y; _A[0]=2;
        assemblyMatrixRow(rowlength);
        for (int m(0); m < rowlength; ++m)
            Matrix.sequential_fill(_column[m], _value[m]);
        Matrix.end_of_row();
        RHS[Zeile] = assemblyRHSLoad(f);
        if(Neumann)
        {
            RHS[Zeile] += assemblyRHSNeumann(1, false, h);
            RHS[Zeile] += assemblyRHSNeumann(2, true, h);
            RHS[Zeile] += assemblyRHSNeumann(3, false, h);
        }
    }
}

    //Kante 2:
    _e.resize(2);
    _A.resize(2);
    rowlength = 12;
    for(int i(1); i<_nx-1; i++)
    {
        Zeile++;
if (Zeile >= fron && Zeile <= lron)
{
        if(Dirichlet)
        {
            Matrix.sequential_fill(Zeile, 1.0);
            Matrix.end_of_row();
            RHS[Zeile]= g.eval(getx(Zeile), gety(Zeile), getz(Zeile));
        }
        else
        {
            _e[0]=Zeile-y-1; _A[0]=3;
            _e[1]=Zeile-y; _A[1]=2;
            assemblyMatrixRow(rowlength);
            for (int m(0); m < rowlength; ++m)
                Matrix.sequential_fill(_column[m], _value[m]);
            Matrix.end_of_row();
            RHS[Zeile] = assemblyRHSLoad(f);
            if(Neumann)
            {
                RHS[Zeile] += assemblyRHSNeumann(1, false, h);
                RHS[Zeile] += assemblyRHSNeumann(2, true, h);
            }
        }
}
    }//close I-Schleife (X-Achse)

    //Ecke 4:
    _e.resize(1);
    _A.resize(1);
    rowlength = 8;
    Zeile++;
if (Zeile >= fron && Zeile <= lron)
{ //Zeile sollte hier z-1 sein
    if(Dirichlet)
    {
        Matrix.sequential_fill(Zeile, 1.0);
        Matrix.end_of_row();
        RHS[Zeile]= g.eval(getx(Zeile), gety(Zeile), getz(Zeile));
    }
    else
    {
        _e[0]=Zeile-y-1; _A[0]=3;
        assemblyMatrixRow(rowlength);
        for (int m(0); m < rowlength; ++m)
            Matrix.sequential_fill(_column[m], _value[m]);
        Matrix.end_of_row();
        RHS[Zeile] = assemblyRHSLoad(f);
        if(Neumann)
        {
            RHS[Zeile] += assemblyRHSNeumann(1, false, h);
            RHS[Zeile] += assemblyRHSNeumann(2, true, h);
            RHS[Zeile] += assemblyRHSNeumann(3, true, h);
        }
    }
}

    for(int k(1); k<_nz-1; k++)
    {
LOG_INFO("assembled ", static_cast<float>(k)/static_cast<double>(_nz)*100.0, "%");
        //Kante 9:
        _e.resize(2);
        _A.resize(2);
        rowlength = 12;
        Zeile++;
if (Zeile >= fron && Zeile <= lron)
{
        if(Dirichlet)
        {
            Matrix.sequential_fill(Zeile, 1.0);
            Matrix.end_of_row();
            RHS[Zeile]= g.eval(getx(Zeile), gety(Zeile), getz(Zeile));
        }
        else
        {
            _e[0]=Zeile-z; _A[0]=4;
            _e[1]=Zeile; _A[1]=0;
            assemblyMatrixRow(rowlength);
            for (int m(0); m < rowlength; ++m)
                Matrix.sequential_fill(_column[m], _value[m]);
            Matrix.end_of_row();
            RHS[Zeile] = assemblyRHSLoad(f);
            if(Neumann)
            {
                RHS[Zeile] += assemblyRHSNeumann(2, false, h);
                RHS[Zeile] += assemblyRHSNeumann(3, false, h);
            }
        }
}

        //Flaeche 3:
        _e.resize(4);
        _A.resize(4);
        rowlength = 18;
        for(int i(1); i<_nx-1; i++)
        {
            Zeile++;
if (Zeile >= fron && Zeile <= lron)
{
            if(Dirichlet)
            {
                Matrix.sequential_fill(Zeile, 1.0);
                Matrix.end_of_row();
                RHS[Zeile]= g.eval(getx(Zeile), gety(Zeile), getz(Zeile));
            }
            else
            {
                _e[0]=Zeile -1-z; _A[0]=5;
                _e[1]=Zeile -z; _A[1]=4;
                _e[2]=Zeile -1; _A[2]=1;
                _e[3]=Zeile; _A[3]=0;
                assemblyMatrixRow(rowlength);
                for (int m(0); m < rowlength; ++m)
                    Matrix.sequential_fill(_column[m], _value[m]);
                Matrix.end_of_row();
                RHS[Zeile] = assemblyRHSLoad(f);
                if(Neumann)
                {
                    RHS[Zeile] += assemblyRHSNeumann(2, false, h);
                }
            }
}
        }//close I-Schleife (X-Achse)

        //Kante 10:
        _e.resize(2);
        _A.resize(2);
        rowlength = 12;
        Zeile++;
if (Zeile >= fron && Zeile <= lron)
{
        if(Dirichlet)
        {
            Matrix.sequential_fill(Zeile, 1.0);
            Matrix.end_of_row();
            RHS[Zeile]= g.eval(getx(Zeile), gety(Zeile), getz(Zeile));
        }
        else
        {
            _e[0]=Zeile-1 -z; _A[0]=5;
            _e[1]=Zeile-1; _A[1]=1;
            assemblyMatrixRow(rowlength);
            for (int m(0); m < rowlength; ++m)
                Matrix.sequential_fill(_column[m], _value[m]);
            Matrix.end_of_row();
            RHS[Zeile] = assemblyRHSLoad(f);
            if(Neumann)
            {
                RHS[Zeile] += assemblyRHSNeumann(2, false, h);
                RHS[Zeile] += assemblyRHSNeumann(3, true, h);
            }
        }
}

        for(int j(1); j<_ny-1; j++)
        {
            //Flaeche 5:
            _e.resize(4);
            _A.resize(4);
            rowlength = 18;
            Zeile++;
if (Zeile >= fron && Zeile <= lron)
{
            if(Dirichlet)
            {
                Matrix.sequential_fill(Zeile, 1.0);
                Matrix.end_of_row();
                RHS[Zeile]= g.eval(getx(Zeile), gety(Zeile), getz(Zeile));
            }
            else
            {
                _e[0]=Zeile -y-z; _A[0]=6;
                _e[1]=Zeile -z; _A[1]=4;
                _e[2]=Zeile -y; _A[2]=2;
                _e[3]=Zeile; _A[3]=0;
                assemblyMatrixRow(rowlength);
                for (int m(0); m < rowlength; ++m)
                    Matrix.sequential_fill(_column[m], _value[m]);
                Matrix.end_of_row();
                RHS[Zeile] = assemblyRHSLoad(f);
                if(Neumann)
                {
                    RHS[Zeile] += assemblyRHSNeumann(3, false, h);
                }
            }
}

            //Inneres:
            for(int i(1); i<_nx-1; i++)
            {
                Zeile++;
if (Zeile >= fron && Zeile <= lron)
{
                if (disc_points[Zeile] == 'o')
                {
                    Matrix.sequential_fill(Zeile, 1.0);
                    Matrix.end_of_row();
                    RHS[Zeile] = -1.0;
                }
                else if (disc_points[Zeile] == 'a')
                {
                    _e.resize(8);
                    _A.resize(8);
                    _e[0]=Zeile -1-y-z; _A[0]=7;
                    _e[1]=Zeile -y-z; _A[1]=6;
                    _e[2]=Zeile -1-z; _A[2]=5;
                    _e[3]=Zeile -z; _A[3]=4;
                    _e[4]=Zeile -1-y; _A[4]=3;
                    _e[5]=Zeile -y; _A[5]=2;
                    _e[6]=Zeile -1; _A[6]=1;
                    _e[7]=Zeile; _A[7]=0;
                    rowlength = 27;

                    assemblyMatrixRow(rowlength);
                    for (int m(0); m < rowlength; ++m)
                        Matrix.sequential_fill(_column[m], _value[m]);
                    Matrix.end_of_row();
                    RHS[Zeile] = assemblyRHSLoad(f);
                }
                else if (disc_points[Zeile] == 'b')
                {
                    if (Dirichlet)
                    {
                        Matrix.sequential_fill(Zeile, 1.0);
                        Matrix.end_of_row();
                        RHS[Zeile]= g.eval(getx(Zeile), gety(Zeile), getz(Zeile));
                    }
                    else
                    {
                        rowlength = setup_A(Zeile, disc_points);
                        setup_e(Zeile);
                        assemblyMatrixRow(rowlength);
                        for (int m(0); m < rowlength; ++m)
                            Matrix.sequential_fill(_column[m], _value[m]);
                        Matrix.end_of_row();
                        RHS[Zeile] = assemblyRHSLoad(f);
                        if (Neumann)
                        {
//                            std::vector<int> planes;
//                            std::vector<bool> rightbacktops;
//                            setup_plane_rigthbacktop(Zeile, disc_points, planes, rightbacktops);
//                            for (int n(0); n < planes.size(); ++n)
//                            {
//                                setup_neumann(Zeile, planes[n], rightbacktops[n], disc_points);
//                                RHS[Zeile] += assemblyRHSNeumann(planes[n], rightbacktops[n], h);
//                            }
                            setup_neumann(Zeile, 1, false, disc_points);
                            RHS[Zeile] += assemblyRHSNeumann(1, false, h);
                            setup_neumann(Zeile, 1, true, disc_points);
                            RHS[Zeile] += assemblyRHSNeumann(1, true, h);
                            setup_neumann(Zeile, 2, false, disc_points);
                            RHS[Zeile] += assemblyRHSNeumann(2, false, h);
                            setup_neumann(Zeile, 2, true, disc_points);
                            RHS[Zeile] += assemblyRHSNeumann(2, true, h);
                            setup_neumann(Zeile, 3, false, disc_points);
                            RHS[Zeile] += assemblyRHSNeumann(3, false, h);
                            setup_neumann(Zeile, 3, true, disc_points);
                            RHS[Zeile] += assemblyRHSNeumann(3, true, h);
                        }
                    }
                }
                else
                {
                    assert(false);
                }
}
            } //Close I-Schleife (X-Achse)

            //Flaeche 6:
            _e.resize(4);
            _A.resize(4);
            rowlength = 18;
            Zeile++;
if (Zeile >= fron && Zeile <= lron)
{
            if(Dirichlet)
            {
                Matrix.sequential_fill(Zeile, 1.0);
                Matrix.end_of_row();
                RHS[Zeile]= g.eval(getx(Zeile), gety(Zeile), getz(Zeile));
            }
            else
            {
                _e[0]=Zeile -y-z -1; _A[0]=7;
                _e[1]=Zeile -z -1; _A[1]=5;
                _e[2]=Zeile -y -1; _A[2]=3;
                _e[3]=Zeile -1; _A[3]=1;
                assemblyMatrixRow(rowlength);
                for (int m(0); m < rowlength; ++m)
                    Matrix.sequential_fill(_column[m], _value[m]);
                Matrix.end_of_row();
                RHS[Zeile] = assemblyRHSLoad(f);
                if(Neumann)
                {
                    RHS[Zeile] += assemblyRHSNeumann(3, true, h);
                }
            }
}
        } //close J-Schleife (Y-Achse)

        //Kante 12:
        _e.resize(2);
        _A.resize(2);
        rowlength = 12;
        Zeile++;
if (Zeile >= fron && Zeile <= lron)
{
        if(Dirichlet)
        {
            Matrix.sequential_fill(Zeile, 1.0);
            Matrix.end_of_row();
            RHS[Zeile]= g.eval(getx(Zeile), gety(Zeile), getz(Zeile));
        }
        else
        {
            _e[0]=Zeile-y -z; _A[0]=6;
            _e[1]=Zeile-y; _A[1]=2;
            assemblyMatrixRow(rowlength);
            for (int m(0); m < rowlength; ++m)
                Matrix.sequential_fill(_column[m], _value[m]);
            Matrix.end_of_row();
            RHS[Zeile] = assemblyRHSLoad(f);
            if(Neumann)
            {
                RHS[Zeile] += assemblyRHSNeumann(2, true, h);
                RHS[Zeile] += assemblyRHSNeumann(3, false, h);
            }
        }
}

        //Flaeche 4
        _e.resize(4);
        _A.resize(4);
        rowlength = 18;
        for(int i(1); i< _nx-1; i++)
        {
            Zeile++;
if (Zeile >= fron && Zeile <= lron)
{
            if(Dirichlet)
            {
                Matrix.sequential_fill( Zeile, 1.0);
                Matrix.end_of_row();
                RHS[Zeile]= g.eval(getx(Zeile), gety(Zeile), getz(Zeile));
            }
            else
            {
                _e[0]= Zeile -1-z -y; _A[0]=7;
                _e[1]= Zeile -z -y; _A[1]=6;
                _e[2]= Zeile -1 -y; _A[2]=3;
                _e[3]= Zeile -y; _A[3]=2;
                assemblyMatrixRow(rowlength);
                for (int m(0); m < rowlength; ++m)
                    Matrix.sequential_fill(_column[m], _value[m]);
                Matrix.end_of_row();
                RHS[Zeile] = assemblyRHSLoad(f);
                if(Neumann)
                {
                    RHS[Zeile] += assemblyRHSNeumann(2, true, h);
                }
            }
}
        }//Close I-Schleife (X-Achse)

        //Kante 11:
        _e.resize(2);
        _A.resize(2);
        rowlength = 12;
        Zeile++;
if (Zeile >= fron && Zeile <= lron)
{
        if(Dirichlet)
        {
            Matrix.sequential_fill(Zeile, 1.0);
            Matrix.end_of_row();
            RHS[Zeile]= g.eval(getx(Zeile), gety(Zeile), getz(Zeile));
        }
        else
        {
            _e[0]=Zeile-1-y-z; _A[0]=7;
            _e[1]=Zeile-1-y; _A[1]=3;
            assemblyMatrixRow(rowlength);
            for (int m(0); m < rowlength; ++m)
                Matrix.sequential_fill(_column[m], _value[m]);
            Matrix.end_of_row();
            RHS[Zeile] = assemblyRHSLoad(f);
            if(Neumann)
            {
                RHS[Zeile] += assemblyRHSNeumann(2, true, h);
                RHS[Zeile] += assemblyRHSNeumann(3, true, h);
            }
        }
}
    } //close K-schleife (Z-Achse)
LOG_INFO("assembled ", static_cast<float>(_nz-1)/static_cast<double>(_nz)*100.0, "%");

    //Ecke 5
    _e.resize(1);
    _A.resize(1);
    rowlength = 8;
    Zeile++;
if (Zeile >= fron && Zeile <= lron)
{ //Zeile sollte hier (_nz-1)*z sein
    if(Dirichlet)
    {
        Matrix.sequential_fill( Zeile, 1.0);
        Matrix.end_of_row();
        RHS[Zeile]= g.eval(getx(Zeile), gety(Zeile), getz(Zeile));
    }
    else
    {
        _e[0]= Zeile-z; _A[0]=4;
        assemblyMatrixRow(rowlength);
        for (int m(0); m < rowlength; ++m)
            Matrix.sequential_fill(_column[m], _value[m]);
        Matrix.end_of_row();
        RHS[Zeile] = assemblyRHSLoad(f);
        if(Neumann)
        {
            RHS[Zeile] += assemblyRHSNeumann(1, true, h);
            RHS[Zeile] += assemblyRHSNeumann(2, false, h);
            RHS[Zeile] += assemblyRHSNeumann(3, false, h);
        }
    }
}

    //Kante 4
    _e.resize(2);
    _A.resize(2);
    rowlength = 12;
    for(int i(1); i<_nx-1;i++)
    {
        Zeile++;
if (Zeile >= fron && Zeile <= lron)
{
        if(Dirichlet)
        {
            Matrix.sequential_fill(Zeile, 1.0);
            Matrix.end_of_row();
            RHS[Zeile]= g.eval(getx(Zeile), gety(Zeile), getz(Zeile));
        }
        else
        {
            _e[0]=Zeile-z -1; _A[0]=5;
            _e[1]=Zeile-z; _A[1]=4;
            assemblyMatrixRow(rowlength);
            for (int m(0); m < rowlength; ++m)
                Matrix.sequential_fill(_column[m], _value[m]);
            Matrix.end_of_row();
            RHS[Zeile] = assemblyRHSLoad(f);
            if(Neumann)
            {
                RHS[Zeile] += assemblyRHSNeumann(1, true, h);
                RHS[Zeile] += assemblyRHSNeumann(2, false, h);
            }
        }
}
    }//Close I-Schleife (X-Achse)

    //Ecke 6:
    _e.resize(1);
    _A.resize(1);
    rowlength = 8;
    Zeile++;
if (Zeile >= fron && Zeile <= lron)
{ //Zeile sollte hier (_nz-1)*z+y-1 sein
    if(Dirichlet)
    {
        Matrix.sequential_fill( Zeile, 1.0);
        Matrix.end_of_row();
        RHS[Zeile]= g.eval(getx(Zeile), gety(Zeile), getz(Zeile));
    }
    else
    {
        _e[0]= Zeile-1-z; _A[0]=5;
        assemblyMatrixRow(rowlength);
        for (int m(0); m < rowlength; ++m)
            Matrix.sequential_fill(_column[m], _value[m]);
        Matrix.end_of_row();
        RHS[Zeile] = assemblyRHSLoad(f);
        if(Neumann)
        {
            RHS[Zeile] += assemblyRHSNeumann(1, true, h);
            RHS[Zeile] += assemblyRHSNeumann(2, false, h);
            RHS[Zeile] += assemblyRHSNeumann(3, true, h);
        }
    }
}

    for(int j(1); j< _ny-1; j++)
    {
        //Kante 8
        _e.resize(2);
        _A.resize(2);
        rowlength = 12;
        Zeile++;
if (Zeile >= fron && Zeile <= lron)
{
        if(Dirichlet)
        {
            Matrix.sequential_fill(Zeile, 1.0);
            Matrix.end_of_row();
            RHS[Zeile]= g.eval(getx(Zeile), gety(Zeile), getz(Zeile));
        }
        else
        {
            _e[0]=Zeile -z-y; _A[0]=6;
            _e[1]=Zeile -z; _A[1]=4;
            assemblyMatrixRow(rowlength);
            for (int m(0); m < rowlength; ++m)
                Matrix.sequential_fill(_column[m], _value[m]);
            Matrix.end_of_row();
            RHS[Zeile] = assemblyRHSLoad(f);
            if(Neumann)
            {
                RHS[Zeile] += assemblyRHSNeumann(1, true, h);
                RHS[Zeile] += assemblyRHSNeumann(3, false, h);
            }
        }
}

        //Flaeche 2:
        _e.resize(4);
        _A.resize(4);
        rowlength = 18;
        for(int i(1); i<_nx-1; i++)
        {
            Zeile++;
if (Zeile >= fron && Zeile <= lron)
{
            if(Dirichlet)
            {
                Matrix.sequential_fill(Zeile, 1.0);
                Matrix.end_of_row();
                RHS[Zeile]= g.eval(getx(Zeile), gety(Zeile), getz(Zeile));
            }
            else
            {
                _e[0]=Zeile -y-1 -z; _A[0]=7;
                _e[1]=Zeile -y -z; _A[1]=6;
                _e[2]=Zeile -1 -z; _A[2]=5;
                _e[3]=Zeile -z ; _A[3]=4;
                assemblyMatrixRow(rowlength);
                for (int m(0); m < rowlength; ++m)
                    Matrix.sequential_fill(_column[m], _value[m]);
                Matrix.end_of_row();
                RHS[Zeile] = assemblyRHSLoad(f);
                if(Neumann)
                {
                    RHS[Zeile] += assemblyRHSNeumann(1, true, h);
                }
            }
}
        }//Close I-Schleife (X-Achse)

        //Kante 7
        _e.resize(2);
        _A.resize(2);
        rowlength = 12;
        Zeile++;
if (Zeile >= fron && Zeile <= lron)
{
        if(Dirichlet)
        {
            Matrix.sequential_fill(Zeile, 1.0);
            Matrix.end_of_row();
            RHS[Zeile]= g.eval(getx(Zeile), gety(Zeile), getz(Zeile));
        }
        else
        {
            _e[0]=Zeile-1-z -y; _A[0]=7;
            _e[1]=Zeile-1-z; _A[1]=5;
            assemblyMatrixRow(rowlength);
            for (int m(0); m < rowlength; ++m)
                Matrix.sequential_fill(_column[m], _value[m]);
            Matrix.end_of_row();
            RHS[Zeile] = assemblyRHSLoad(f);
            if(Neumann)
            {
                RHS[Zeile] += assemblyRHSNeumann(1, true, h);
                RHS[Zeile] += assemblyRHSNeumann(3, true, h);
            }
        }
}
    }//Close J-Schleife (Y-Achse)

    //Ecke 7
    _e.resize(1);
    _A.resize(1);
    rowlength = 8;
    Zeile++;
if (Zeile >= fron && Zeile <= lron)
{ //Zeile sollte hier (_nx*_ny*_nz)-_nx sein
    if(Dirichlet)
    {
        Matrix.sequential_fill(Zeile, 1.0);
        Matrix.end_of_row();
        RHS[Zeile]= g.eval(getx(Zeile), gety(Zeile), getz(Zeile));
    }
    else
    {
        _e[0]=Zeile-y-z; _A[0]=6;
        assemblyMatrixRow(rowlength);
        for (int m(0); m < rowlength; ++m)
            Matrix.sequential_fill(_column[m], _value[m]);
        Matrix.end_of_row();
        RHS[Zeile] = assemblyRHSLoad(f);
        if(Neumann)
        {
            RHS[Zeile] += assemblyRHSNeumann(1, true, h);
            RHS[Zeile] += assemblyRHSNeumann(2, true, h);
            RHS[Zeile] += assemblyRHSNeumann(3, false, h);
        }
    }
}

    //Kante 3
    _e.resize(2);
    _A.resize(2);
    rowlength = 12;
    for(int i(1); i<_nx-1; i++)
    {
        Zeile++;
if (Zeile >= fron && Zeile <= lron)
{
        if(Dirichlet)
        {
            Matrix.sequential_fill( Zeile, 1.0);
            Matrix.end_of_row();
            RHS[Zeile]= g.eval(getx(Zeile), gety(Zeile), getz(Zeile));
        }
        else
        {
            _e[0]= Zeile-y-z -1; _A[0]=7;
            _e[1]= Zeile-y-z; _A[1]=6;
            assemblyMatrixRow(rowlength);
            for (int m(0); m < rowlength; ++m)
                Matrix.sequential_fill(_column[m], _value[m]);
            Matrix.end_of_row();
            RHS[Zeile] = assemblyRHSLoad(f);
            if(Neumann)
            {
                RHS[Zeile] += assemblyRHSNeumann(1, true, h);
                RHS[Zeile] += assemblyRHSNeumann(2, true, h);
            }
        }
}
    }//Close I-Schleife (X-Achse)

    //Ecke 8:
    _e.resize(1);
    _A.resize(1);
    rowlength = 8;
    Zeile++;
if (Zeile >= fron && Zeile <= lron)
{ //Zeile sollte hier (_nx*_ny*_nz) sein
    if(Dirichlet)
    {
        Matrix.sequential_fill(Zeile, 1.0);
        Matrix.end_of_row();
        RHS[Zeile]= g.eval(getx(Zeile), gety(Zeile), getz(Zeile));
    }
    else
    {
        _e[0]=Zeile-z-y-1; _A[0]=7;
        assemblyMatrixRow(rowlength);
        for (int m(0); m < rowlength; ++m)
            Matrix.sequential_fill(_column[m], _value[m]);
        Matrix.end_of_row();
        RHS[Zeile] = assemblyRHSLoad(f);
        if(Neumann)
        {
            RHS[Zeile] += assemblyRHSNeumann(1, true, h);
            RHS[Zeile] += assemblyRHSNeumann(2, true, h);
            RHS[Zeile] += assemblyRHSNeumann(3, true, h);
        }
    }
}
LOG_INFO("assembled 100%");

    //TODO rhs direkt fuellen (erst wenn alles laeuft)
    for (int i(0); i<_nx*_ny*_nz; ++i)
        rhs.set_global(i, RHS[i]);

    LOG_INFO("Matrix succesfully assembled");

}//assemble()

}//namespace Icarus
