#include "include/assemblefem.hpp"
#include "include/mathfunction.hpp"

#include "include/logger.hpp"
#include <iostream>

namespace Icarus
{

void assembleFem::assemble(DistEllpackMatrix<double>& Matrix, SlicedVector<double>& rhs,
    std::vector<char>& disc_points, mathfunction f, mathfunction g, mathfunction h)
{
    //TODO: vorlaeufig, wieder loeschen
    bool Dirichlet(true);
    bool Neumann(false);
    //TODO: vorlaeufig, wieder loeschen

    Matrix.prepare_sequential_fill(27);

    int Zeile;
    std::vector<double> RHS(_nx*_ny*_nz);

LOG_INFO("assembled 0%");
    //Ecke 1
    _e.clear(); _e.resize(1);
    _A.clear(); _A.resize(1);
    Zeile=0;
    if(Dirichlet)
    {
        Matrix.sequential_fill(Zeile, 1.0);
        Matrix.end_of_row();
        RHS[Zeile]= g.eval(getx(Zeile), gety(Zeile), getz(Zeile));
    }
    else
    {
        _e[0]=Zeile; _A[0]=0;
        assemblyMatrixRow();
        for (int m(0); m<8; ++m)
            Matrix.sequential_fill(_column[m], _value[m]);
        Matrix.end_of_row();
        RHS[Zeile] = assemblyRHSLoad(f);
        if(Neumann)
        {
            RHS[Zeile] += assemblyRHSNeumann(1, 0, h);
            RHS[Zeile] += assemblyRHSNeumann(2, 0, h);
            RHS[Zeile] += assemblyRHSNeumann(3, 0, h);
        }
    }

    //Kante 1:
    _e.resize(2);
    _A.resize(2);
    for(int i(1); i<_nx-1;i++)
    {
        Zeile++;
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
            assemblyMatrixRow();
            for (int m(0); m<12; ++m)
                Matrix.sequential_fill(_column[m], _value[m]);
            Matrix.end_of_row();
            RHS[Zeile] = assemblyRHSLoad(f);
            if(Neumann)
            {
                RHS[Zeile] += assemblyRHSNeumann(1, 0, h);
                RHS[Zeile] += assemblyRHSNeumann(2, 0, h);
            }
        }
    }//close I-Schleife (X-Achse)

    //Ecke 2
    _e.resize(1);
    _A.resize(1);
    Zeile++; //Zeile sollte hier y-1 sein
    if(Dirichlet)
    {
        Matrix.sequential_fill(Zeile, 1.0);
        Matrix.end_of_row();
        RHS[Zeile]= g.eval(getx(Zeile), gety(Zeile), getz(Zeile));
    }
    else
    {
        _e[0]=Zeile-1; _A[0]=1;
        assemblyMatrixRow();
        for (int m(0); m<8; ++m)
            Matrix.sequential_fill(_column[m], _value[m]);
        Matrix.end_of_row();
        RHS[Zeile] = assemblyRHSLoad(f);
        if(Neumann)
        {
            RHS[Zeile] += assemblyRHSNeumann(1, 0, h);
            RHS[Zeile] += assemblyRHSNeumann(2, 0, h);
            RHS[Zeile] += assemblyRHSNeumann(3, 1, h);
        }
    }

    for(int j(1); j<_ny-1;j++)
    {
        //Kante 5
        _e.resize(2);
        _A.resize(2);

        Zeile++;
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
            assemblyMatrixRow();
            for (int m(0); m<12; ++m)
                Matrix.sequential_fill(_column[m], _value[m]);
            Matrix.end_of_row();
            RHS[Zeile] = assemblyRHSLoad(f);
            if(Neumann)
            {
                RHS[Zeile] += assemblyRHSNeumann(1, 0, h);
                RHS[Zeile] += assemblyRHSNeumann(3, 0, h);
            }
        }

        //Fläche 1
        _e.resize(4);
        _A.resize(4);
        for(int i(1); i<_nx-1;i++)
        {
            Zeile++;
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
                assemblyMatrixRow();
                for (int m(0); m<18; ++m)
                    Matrix.sequential_fill(_column[m], _value[m]);
                Matrix.end_of_row();
                RHS[Zeile] = assemblyRHSLoad(f);
                if(Neumann)
                {
                    RHS[Zeile] += assemblyRHSNeumann(1, 0, h);
                }
            }
        } //close I-Schleife (X-Achse)

        //Kante: 6
        _e.resize(2);
        _A.resize(2);
        Zeile++;
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
            assemblyMatrixRow();
            for (int m(0); m<12; ++m)
                Matrix.sequential_fill(_column[m], _value[m]);
            Matrix.end_of_row();
            RHS[Zeile] = assemblyRHSLoad(f);
            if(Neumann)
            {
                RHS[Zeile] += assemblyRHSNeumann(1, 0, h);
                RHS[Zeile] += assemblyRHSNeumann(3, 1, h);
            }
        }
    } //close J-Schleife (Y-Achse)

    //Ecke 3:
    _e.resize(1);
    _A.resize(1);
    Zeile++; //Zeile sollte hier (_ny-1)*y sein
    if(Dirichlet)
    {
        Matrix.sequential_fill(Zeile, 1.0);
        Matrix.end_of_row();
        RHS[Zeile]= g.eval(getx(Zeile), gety(Zeile), getz(Zeile));
    }
    else
    {
        _e[0]=Zeile-y; _A[0]=2;
        assemblyMatrixRow();
        for (int m(0); m<8; ++m)
            Matrix.sequential_fill(_column[m], _value[m]);
        Matrix.end_of_row();
        RHS[Zeile] = assemblyRHSLoad(f);
        if(Neumann)
        {
            RHS[Zeile] += assemblyRHSNeumann(1, 0, h);
            RHS[Zeile] += assemblyRHSNeumann(2, 1, h);
            RHS[Zeile] += assemblyRHSNeumann(3, 0, h);
        }
    }

    //Kante 2:
    _e.resize(2);
    _A.resize(2);
    for(int i(1); i<_nx-1; i++)
    {
        Zeile++;
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
            assemblyMatrixRow();
            for (int m(0); m<12; ++m)
                Matrix.sequential_fill(_column[m], _value[m]);
            Matrix.end_of_row();
            RHS[Zeile] = assemblyRHSLoad(f);
            if(Neumann)
            {
                RHS[Zeile] += assemblyRHSNeumann(1, 0, h);
                RHS[Zeile] += assemblyRHSNeumann(2, 1, h);
            }
        }
    }//close I-Schleife (X-Achse)

    //Ecke 4:
    _e.resize(1);
    _A.resize(1);
    Zeile++; //Zeile sollte hier z-1 sein
    if(Dirichlet)
    {
        Matrix.sequential_fill(Zeile, 1.0);
        Matrix.end_of_row();
        RHS[Zeile]= g.eval(getx(Zeile), gety(Zeile), getz(Zeile));
    }
    else
    {
        _e[0]=Zeile-y-1; _A[0]=3;
        assemblyMatrixRow();
        for (int m(0); m<8; ++m)
            Matrix.sequential_fill(_column[m], _value[m]);
        Matrix.end_of_row();
        RHS[Zeile] = assemblyRHSLoad(f);
        if(Neumann)
        {
            RHS[Zeile] += assemblyRHSNeumann(1, 0, h);
            RHS[Zeile] += assemblyRHSNeumann(2, 1, h);
            RHS[Zeile] += assemblyRHSNeumann(3, 1, h);
        }
    }

    for(int k(1); k<_nz-1; k++)
    {
LOG_INFO("assembled ", static_cast<float>(k)/static_cast<double>(_nz)*100.0, "%");
        //Kante 9:
        _e.resize(2);
        _A.resize(2);
        Zeile++;
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
            assemblyMatrixRow();
            for (int m(0); m<12; ++m)
                Matrix.sequential_fill(_column[m], _value[m]);
            Matrix.end_of_row();
            RHS[Zeile] = assemblyRHSLoad(f);
            if(Neumann)
            {
                RHS[Zeile] += assemblyRHSNeumann(2, 0, h);
                RHS[Zeile] += assemblyRHSNeumann(3, 0, h);
            }
        }

        //Flaeche 3:
        _e.resize(4);
        _A.resize(4);
        for(int i(1); i<_nx-1; i++)
        {
            Zeile++;
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
                assemblyMatrixRow();
                for (int m(0); m<18; ++m)
                    Matrix.sequential_fill(_column[m], _value[m]);
                Matrix.end_of_row();
                RHS[Zeile] = assemblyRHSLoad(f);
                if(Neumann)
                {
                    RHS[Zeile] += assemblyRHSNeumann(2, 0, h);
                }
            }
        }//close I-Schleife (X-Achse)

        //Kante 10:
        _e.resize(2);
        _A.resize(2);
        Zeile++;
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
            assemblyMatrixRow();
            for (int m(0); m<12; ++m)
                Matrix.sequential_fill(_column[m], _value[m]);
            Matrix.end_of_row();
            RHS[Zeile] = assemblyRHSLoad(f);
            if(Neumann)
            {
                RHS[Zeile] += assemblyRHSNeumann(2, 0, h);
                RHS[Zeile] += assemblyRHSNeumann(3, 1, h);
            }
        }

        for(int j(1); j<_ny-1; j++)
        {
            //Fläche 5:
            _e.resize(4);
            _A.resize(4);
            Zeile++;
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
                assemblyMatrixRow();
                for (int m(0); m<18; ++m)
                    Matrix.sequential_fill(_column[m], _value[m]);
                Matrix.end_of_row();
                RHS[Zeile] = assemblyRHSLoad(f);
                if(Neumann)
                {
                    RHS[Zeile] += assemblyRHSNeumann(3, 0, h);
                }
            }

            //Inneres:
            _e.resize(8);
            _A.resize(8);
            for(int i(1); i<_nx-1; i++)
            {
                Zeile++;;
                //if(Dirichlet)
                //{
                //    Matrix.sequential_fill(Zeile, 1.0);
                //    Matrix.end_of_row();
                //    RHS[Zeile]= g.eval(getx(Zeile), gety(Zeile), getz(Zeile));
                //}
                //else
                {
                    _e[0]=Zeile -1-y-z; _A[0]=7;
                    _e[1]=Zeile -y-z; _A[1]=6;
                    _e[2]=Zeile -1-z; _A[2]=5;
                    _e[3]=Zeile -z; _A[3]=4;
                    _e[4]=Zeile -1-y; _A[4]=3;
                    _e[5]=Zeile -y; _A[5]=2;
                    _e[6]=Zeile -1; _A[6]=1;
                    _e[7]=Zeile; _A[7]=0;

                    assemblyMatrixRow();
                    for (int m(0); m<27; ++m)
                        Matrix.sequential_fill(_column[m], _value[m]);
                    Matrix.end_of_row();
                    RHS[Zeile] = assemblyRHSLoad(f);

                    //Neumann eventuell hinzufuegen
                }
            } //Close I-Schleife (X-Achse)

            //Fläche 6:
            _e.resize(4);
            _A.resize(4);
            Zeile++;
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
                assemblyMatrixRow();
                for (int m(0); m<18; ++m)
                    Matrix.sequential_fill(_column[m], _value[m]);
                Matrix.end_of_row();
                RHS[Zeile] = assemblyRHSLoad(f);
                if(Neumann)
                {
                    RHS[Zeile] += assemblyRHSNeumann(3, 1, h);
                }
            }
        } //close J-Schleife (Y-Achse)

        //Kante 12:
        _e.resize(2);
        _A.resize(2);
        Zeile++;
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
            assemblyMatrixRow();
            for (int m(0); m<12; ++m)
                Matrix.sequential_fill(_column[m], _value[m]);
            Matrix.end_of_row();
            RHS[Zeile] = assemblyRHSLoad(f);
            if(Neumann)
            {
                RHS[Zeile] += assemblyRHSNeumann(2, 1, h);
                RHS[Zeile] += assemblyRHSNeumann(3, 0, h);
            }
        }

        //Fläche 4
        _e.resize(4);
        _A.resize(4);
        for(int i(1); i< _nx-1; i++)
        {
            Zeile++;
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
                assemblyMatrixRow();
                for (int m(0); m<18; ++m)
                    Matrix.sequential_fill(_column[m], _value[m]);
                Matrix.end_of_row();
                RHS[Zeile] = assemblyRHSLoad(f);
                if(Neumann)
                {
                    RHS[Zeile] += assemblyRHSNeumann(2, 1, h);
                }
            }
        }//Close I-Schleife (X-Achse)

        //Kante 11:
        _e.resize(2);
        _A.resize(2);
        Zeile++;
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
            assemblyMatrixRow();
            for (int m(0); m<12; ++m)
                Matrix.sequential_fill(_column[m], _value[m]);
            Matrix.end_of_row();
            RHS[Zeile] = assemblyRHSLoad(f);
            if(Neumann)
            {
                RHS[Zeile] += assemblyRHSNeumann(2, 1, h);
                RHS[Zeile] += assemblyRHSNeumann(3, 1, h);
            }
        }
    } //close K-schleife (Z-Achse)
LOG_INFO("assembled ", static_cast<float>(_nz-1)/static_cast<double>(_nz)*100.0, "%");

    //Ecke 5
    _e.resize(1);
    _A.resize(1);
    Zeile++; //Zeile sollte hier (_nz-1)*z sein
    if(Dirichlet)
    {
        Matrix.sequential_fill( Zeile, 1.0);
        Matrix.end_of_row();
        RHS[Zeile]= g.eval(getx(Zeile), gety(Zeile), getz(Zeile));
    }
    else
    {
        _e[0]= Zeile-z; _A[0]=4;
        assemblyMatrixRow();
        for (int m(0); m<8; ++m)
            Matrix.sequential_fill(_column[m], _value[m]);
        Matrix.end_of_row();
        RHS[Zeile] = assemblyRHSLoad(f);
        if(Neumann)
        {
            RHS[Zeile] += assemblyRHSNeumann(1, 1, h);
            RHS[Zeile] += assemblyRHSNeumann(2, 0, h);
            RHS[Zeile] += assemblyRHSNeumann(3, 0, h);
        }
    }

    //Kante 4
    _e.resize(2);
    _A.resize(2);
    for(int i(1); i<_nx-1;i++)
    {
        Zeile++;
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
            assemblyMatrixRow();
            for (int m(0); m<12; ++m)
                Matrix.sequential_fill(_column[m], _value[m]);
            Matrix.end_of_row();
            RHS[Zeile] = assemblyRHSLoad(f);
            if(Neumann)
            {
                RHS[Zeile] += assemblyRHSNeumann(1, 1, h);
                RHS[Zeile] += assemblyRHSNeumann(2, 0, h);
            }
        }
    }//Close I-Schleife (X-Achse)

    //Ecke 6:
    _e.resize(1);
    _A.resize(1);
    Zeile++; //Zeile sollte hier (_nz-1)*z+y-1 sein
    if(Dirichlet)
    {
        Matrix.sequential_fill( Zeile, 1.0);
        Matrix.end_of_row();
        RHS[Zeile]= g.eval(getx(Zeile), gety(Zeile), getz(Zeile));
    }
    else
    {
        _e[0]= Zeile-1-z; _A[0]=5;
        assemblyMatrixRow();
        for (int m(0); m<8; ++m)
            Matrix.sequential_fill(_column[m], _value[m]);
        Matrix.end_of_row();
        RHS[Zeile] = assemblyRHSLoad(f);
        if(Neumann)
        {
            RHS[Zeile] += assemblyRHSNeumann(1, 1, h);
            RHS[Zeile] += assemblyRHSNeumann(2, 0, h);
            RHS[Zeile] += assemblyRHSNeumann(3, 1, h);
        }
    }

    for(int j(1); j< _ny-1; j++)
    {
        //Kante 8
        _e.resize(2);
        _A.resize(2);
        Zeile++;
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
            assemblyMatrixRow();
            for (int m(0); m<12; ++m)
                Matrix.sequential_fill(_column[m], _value[m]);
            Matrix.end_of_row();
            RHS[Zeile] = assemblyRHSLoad(f);
            if(Neumann)
            {
                RHS[Zeile] += assemblyRHSNeumann(1, 1, h);
                RHS[Zeile] += assemblyRHSNeumann(3, 0, h);
            }
        }

        //Fläche 2:
        _e.resize(4);
        _A.resize(4);
        for(int i(1); i<_nx-1; i++)
        {
            Zeile++;
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
                assemblyMatrixRow();
                for (int m(0); m<18; ++m)
                    Matrix.sequential_fill(_column[m], _value[m]);
                Matrix.end_of_row();
                RHS[Zeile] = assemblyRHSLoad(f);
                if(Neumann)
                {
                    RHS[Zeile] += assemblyRHSNeumann(1, 1, h);
                }
            }
        }//Close I-Schleife (X-Achse)

        //Kante 7
        _e.resize(2);
        _A.resize(2);
        Zeile++;
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
            assemblyMatrixRow();
            for (int m(0); m<12; ++m)
                Matrix.sequential_fill(_column[m], _value[m]);
            Matrix.end_of_row();
            RHS[Zeile] = assemblyRHSLoad(f);
            if(Neumann)
            {
                RHS[Zeile] += assemblyRHSNeumann(1, 1, h);
                RHS[Zeile] += assemblyRHSNeumann(3, 1, h);
            }
        }
    }//Close J-Schleife (Y-Achse)

    //Ecke 7
    _e.resize(1);
    _A.resize(1);
    Zeile++; //Zeile sollte hier (_nx*_ny*_nz)-_nx sein
    if(Dirichlet)
    {
        Matrix.sequential_fill(Zeile, 1.0);
        Matrix.end_of_row();
        RHS[Zeile]= g.eval(getx(Zeile), gety(Zeile), getz(Zeile));
    }
    else
    {
        _e[0]=Zeile-y-z; _A[0]=6;
        assemblyMatrixRow();
        for (int m(0); m<8; ++m)
            Matrix.sequential_fill(_column[m], _value[m]);
        Matrix.end_of_row();
        RHS[Zeile] = assemblyRHSLoad(f);
        if(Neumann)
        {
            RHS[Zeile] += assemblyRHSNeumann(1, 1, h);
            RHS[Zeile] += assemblyRHSNeumann(2, 1, h);
            RHS[Zeile] += assemblyRHSNeumann(3, 0, h);
        }
    }

    //Kante 3
    _e.resize(2);
    _A.resize(2);
    for(int i(1); i<_nx-1; i++)
    {
        Zeile++;
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
            assemblyMatrixRow();
            for (int m(0); m<12; ++m)
                Matrix.sequential_fill(_column[m], _value[m]);
            Matrix.end_of_row();
            RHS[Zeile] = assemblyRHSLoad(f);
            if(Neumann)
            {
                RHS[Zeile] += assemblyRHSNeumann(1, 1, h);
                RHS[Zeile] += assemblyRHSNeumann(2, 1, h);
            }
        }
    }//Close I-Schleife (X-Achse)

    //Ecke 8:
    _e.resize(1);
    _A.resize(1);
    Zeile++; //Zeile sollte hier (_nx*_ny*_nz) sein
    if(Dirichlet)
    {
        Matrix.sequential_fill(Zeile, 1.0);
        Matrix.end_of_row();
        RHS[Zeile]= g.eval(getx(Zeile), gety(Zeile), getz(Zeile));
    }
    else
    {
        _e[0]=Zeile-z-y-1; _A[0]=7;
        assemblyMatrixRow();
        for (int m(0); m<8; ++m)
            Matrix.sequential_fill(_column[m], _value[m]);
        Matrix.end_of_row();
        RHS[Zeile] = assemblyRHSLoad(f);
        if(Neumann)
        {
            RHS[Zeile] += assemblyRHSNeumann(1, 1, h);
            RHS[Zeile] += assemblyRHSNeumann(2, 1, h);
            RHS[Zeile] += assemblyRHSNeumann(3, 1, h);
        }
    }
LOG_INFO("assembled 100%");

    //TODO rhs direkt fuellen (erst wenn alles laeuft)
    for (int i(0); i<_nx*_ny*_nz; ++i)
        rhs.set_global(i, RHS[i]);

    LOG_INFO("Matrix succesfully assembled");

}//assemble()

}//namespace Icarus
