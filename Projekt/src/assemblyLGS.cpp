#include "include/assemblefem.hpp"
#include "include/mathfunction.hpp"

#include <iostream>

namespace Icarus
{

void assembleFem::assemble(DistEllpackMatrix<double>& Matrix, SlicedVector<double>& rhs, mathfunction f)
{
    //TODO: vorlaeufig, wieder loeschen
    bool Dirichlet(true);
    bool Neumann(false);
    double RHSVAL(1.0);
    //TODO: vorlaeufig, wieder loeschen

    Matrix.prepare_sequential_fill(27);

    int Zeile;
    std::vector<int> e(1);
    std::vector<int> A(1);
    std::vector<double> RHS(_nx*_ny*_nz);

    std::vector<int> column(27);
    std::vector<double> value(27);

    //Ecke 1
    Zeile=0;
    if(Dirichlet)
    {
        Matrix.sequential_fill(Zeile, 1.0);
        Matrix.end_of_row();
        RHS[Zeile]= RHSVAL;
    }
    else
    {
        e[0]=Zeile; A[0]=0;
        assemblyMatrixRow(e, A, column, value);
        for (int m(0); m<8; ++m)
            Matrix.sequential_fill(column[m], value[m]);
        Matrix.end_of_row();
        RHS[Zeile] = assemblyRHSLoad(e, A, f);
        if(Neumann)
        {
            e[0]=Zeile; A[0]=0;
            RHS[Zeile] += assemblyRHSNeumann(e, A, 1);
            e[0]=Zeile; A[0]=0;
            RHS[Zeile] += assemblyRHSNeumann(e, A, 2);
            e[0]=Zeile; A[0]=0;
            RHS[Zeile] += assemblyRHSNeumann(e, A, 3);
        }
    }

    //Kante 1:
    e.resize(2);
    A.resize(2);

    for(int i=1; i<_nx-1;i++)
    {
        Zeile++;;
        if(Dirichlet)
        {
            Matrix.sequential_fill(Zeile, 1.0);
            Matrix.end_of_row();
            RHS[Zeile]= RHSVAL;
        }
        else
        {
            e[0]=Zeile-1; A[0]=1;
            e[1]=Zeile; A[1]=0;
            assemblyMatrixRow(e, A, column, value);
            for (int m(0); m<12; ++m)
                Matrix.sequential_fill(column[m], value[m]);
            Matrix.end_of_row();
            RHS[Zeile] = assemblyRHSLoad(e, A, f);
            if(Neumann)
            {
                e[0]=Zeile-1; A[0]=1;
                e[1]=Zeile; A[1]=0;
                RHS[Zeile] += assemblyRHSNeumann(e, A, 1);
                e[0]=Zeile-1; A[0]=1;
                e[1]=Zeile; A[1]=0;
                RHS[Zeile] += assemblyRHSNeumann(e, A, 2);
            }
        }
    }//close I-Schleife (X-Achse)

    //Ecke 2
    e.resize(1);
    A.resize(1);
    Zeile++; //Zeile sollte hier y-1 sein
    if(Dirichlet)
    {
        Matrix.sequential_fill(Zeile, 1.0);
        Matrix.end_of_row();
        RHS[Zeile]= RHSVAL;
    }
    else
    {
        e[0]=Zeile-1; A[0]=1;
        assemblyMatrixRow(e, A, column, value);
        for (int m(0); m<8; ++m)
            Matrix.sequential_fill(column[m], value[m]);
        Matrix.end_of_row();
        RHS[Zeile] = assemblyRHSLoad(e, A, f);
        if(Neumann)
        {
            e[0]=Zeile-1; A[0]=1;
            RHS[Zeile] += assemblyRHSNeumann(e, A, 1);
            e[0]=Zeile-1; A[0]=1;
            RHS[Zeile] += assemblyRHSNeumann(e, A, 2);
            e[0]=Zeile; A[0]=0;
            RHS[Zeile] += assemblyRHSNeumann(e, A, 3);
        }
    }

    for(int j=1; j<_ny-1;j++)
    {
        //Kante 5
        e.resize(2);
        A.resize(2);

        Zeile++;
        if(Dirichlet)
        {
            Matrix.sequential_fill(Zeile, 1.0);
            Matrix.end_of_row();
            RHS[Zeile]= RHSVAL;
        }
        else
        {
            e[0]=Zeile-y; A[0]=3;
            e[1]=Zeile; A[1]=0;
            assemblyMatrixRow(e, A, column, value);
            for (int m(0); m<12; ++m)
                Matrix.sequential_fill(column[m], value[m]);
            Matrix.end_of_row();
            RHS[Zeile] = assemblyRHSLoad(e, A, f);
            if(Neumann)
            {
                e[0]=Zeile-y; A[0]=3;
                e[1]=Zeile; A[1]=0;
                RHS[Zeile] += assemblyRHSNeumann(e, A, 1);
                e[0]=Zeile-y; A[0]=1;
                e[1]=Zeile; A[1]=0;
                RHS[Zeile] += assemblyRHSNeumann(e, A, 3);
            }
        }

        //Fläche 1:
        e.resize(4);
        A.resize(4);
        for(int i=1; i<_nx-1;i++)
        {
            Zeile++;
            if(Dirichlet)
            {
                Matrix.sequential_fill(Zeile, 1.0);
                Matrix.end_of_row();
                RHS[Zeile]= RHSVAL;
            }
            else
            {
                e[0]=Zeile -y-1; A[0]=2;
                e[1]=Zeile -y; A[1]=3;
                e[2]=Zeile; A[2]=0;
                e[3]=Zeile -1; A[3]=1;
                assemblyMatrixRow(e, A, column, value);
                for (int m(0); m<18; ++m)
                    Matrix.sequential_fill(column[m], value[m]);
                Matrix.end_of_row();
                RHS[Zeile] = assemblyRHSLoad(e, A, f);
                if(Neumann)
                {
                    e[0]=Zeile -y-1; A[0]=2;
                    e[1]=Zeile -y; A[1]=3;
                    e[2]=Zeile; A[2]=0;
                    e[3]=Zeile -1; A[3]=1;
                    RHS[Zeile] += assemblyRHSNeumann(e, A, 1);
                }
            }
        } //close I-Schleife (X-Achse)

        //Kante: 6
        e.resize(2);
        A.resize(2);
        Zeile++;
        if(Dirichlet)
        {
            Matrix.sequential_fill(Zeile, 1.0);
            Matrix.end_of_row();
            RHS[Zeile]= RHSVAL;
        }
        else
        {
            e[0]=Zeile -1 -y; A[0]=2;
            e[1]=Zeile -1; A[1]=1;
            assemblyMatrixRow(e, A, column, value);
            for (int m(0); m<12; ++m)
                Matrix.sequential_fill(column[m], value[m]);
            Matrix.end_of_row();
            RHS[Zeile] = assemblyRHSLoad(e, A, f);
            if(Neumann)
            {
                e[0]=Zeile-y-1; A[0]=2;
                e[1]=Zeile-1; A[1]=1;
                RHS[Zeile] += assemblyRHSNeumann(e, A, 1);
                e[0]=Zeile-y; A[0]=1;
                e[1]=Zeile; A[1]=0;
                RHS[Zeile] += assemblyRHSNeumann(e, A, 3);
            }
        }
    } //close J-Schleife (Y-Achse)

    //Ecke 3:
    e.resize(1);
    A.resize(1);
    Zeile++; //Zeile sollte hier (_ny-1)*y sein
    if(Dirichlet)
    {
        Matrix.sequential_fill(Zeile, 1.0);
        Matrix.end_of_row();
        RHS[Zeile]= RHSVAL;
    }
    else
    {
        e[0]=Zeile-y; A[0]=3;
        assemblyMatrixRow(e, A, column, value);
        for (int m(0); m<8; ++m)
            Matrix.sequential_fill(column[m], value[m]);
        Matrix.end_of_row();
        RHS[Zeile] = assemblyRHSLoad(e, A, f);
        if(Neumann)
        {
            e[0]=Zeile-y; A[0]=3;
            RHS[Zeile] += assemblyRHSNeumann(e, A, 1);
            e[0]=Zeile; A[0]=0;
            RHS[Zeile] += assemblyRHSNeumann(e, A, 2);
            e[0]=Zeile-y; A[0]=1;
            RHS[Zeile] += assemblyRHSNeumann(e, A, 3);
        }
    }

    //Kante 2:
    e.resize(2);
    A.resize(2);
    for(int i=1; i<_nx-1; i++)
    {
        Zeile++;
        if(Dirichlet)
        {
            Matrix.sequential_fill(Zeile, 1.0);
            Matrix.end_of_row();
            RHS[Zeile]= RHSVAL;
        }
        else
        {
            e[0]=Zeile-y -1; A[0]=2;
            e[1]=Zeile-y; A[1]=3;
            assemblyMatrixRow(e, A, column, value);
            for (int m(0); m<12; ++m)
                Matrix.sequential_fill(column[m], value[m]);
            Matrix.end_of_row();
            RHS[Zeile] = assemblyRHSLoad(e, A, f);
            if(Neumann)
            {
                e[0]=Zeile-y-1; A[0]=2;
                e[1]=Zeile-y; A[1]=3;
                RHS[Zeile] += assemblyRHSNeumann(e, A, 1);
                e[0]=Zeile-1; A[0]=1;
                e[1]=Zeile; A[1]=0;
                RHS[Zeile] += assemblyRHSNeumann(e, A, 2);
            }
        }
    }//close I-Schleife (X-Achse)

    //Ecke 4:
    e.resize(1);
    A.resize(1);
    Zeile++; //Zeile sollte hier z-1 sein
    if(Dirichlet)
    {
        Matrix.sequential_fill(Zeile, 1.0);
        Matrix.end_of_row();
        RHS[Zeile]= RHSVAL;
    }
    else
    {
        e[0]=Zeile-1-y; A[0]=2;
        assemblyMatrixRow(e, A, column, value);
        for (int m(0); m<8; ++m)
            Matrix.sequential_fill(column[m], value[m]);
        Matrix.end_of_row();
        RHS[Zeile] = assemblyRHSLoad(e, A, f);
        if(Neumann)
        {
            e[0]=Zeile-1-y; A[0]=2;
            RHS[Zeile] += assemblyRHSNeumann(e, A, 1);
            e[0]=Zeile-1; A[0]=1;
            RHS[Zeile] += assemblyRHSNeumann(e, A, 2);
            e[0]=Zeile-y; A[0]=1;
            RHS[Zeile] += assemblyRHSNeumann(e, A, 3);
        }
    }

    for(int k=1; k<_nz-1; k++)
    {
        //Kante 9:
        e.resize(2);
        A.resize(2);
        Zeile++;
        if(Dirichlet)
        {
            Matrix.sequential_fill(Zeile, 1.0);
            Matrix.end_of_row();
            RHS[Zeile]= RHSVAL;
        }
        else
        {
            e[0]=Zeile-z; A[0]=4;
            e[1]=Zeile; A[1]=0;
            assemblyMatrixRow(e, A, column, value);
            for (int m(0); m<12; ++m)
                Matrix.sequential_fill(column[m], value[m]);
            Matrix.end_of_row();
            RHS[Zeile] = assemblyRHSLoad(e, A, f);
            if(Neumann)
            {
                e[0]=Zeile-z; A[0]=3;
                e[1]=Zeile; A[1]=0;
                RHS[Zeile] += assemblyRHSNeumann(e, A, 2);
                e[0]=Zeile-z; A[0]=3;
                e[1]=Zeile; A[1]=0;
                RHS[Zeile] += assemblyRHSNeumann(e, A, 3);
            }
        }

        //Flaeche 3:
        e.resize(4);
        A.resize(4);
        for(int i=1; i<_nx-1; i++)
        {
            Zeile++;
            if(Dirichlet)
            {
                Matrix.sequential_fill(Zeile, 1.0);
                Matrix.end_of_row();
                RHS[Zeile]= RHSVAL;
            }
            else
            {
                e[0]=Zeile -1-z; A[0]=5;
                e[1]=Zeile -z; A[1]=4;
                e[2]=Zeile; A[2]=0;
                e[3]=Zeile -1; A[3]=1;
                assemblyMatrixRow(e, A, column, value);
                for (int m(0); m<18; ++m)
                    Matrix.sequential_fill(column[m], value[m]);
                Matrix.end_of_row();
                RHS[Zeile] = assemblyRHSLoad(e, A, f);
                if(Neumann)
                {
                    e[0]=Zeile -1-z; A[0]=2;
                    e[1]=Zeile -z; A[1]=3;
                    e[2]=Zeile; A[2]=0;
                    e[3]=Zeile -1; A[3]=1;
                    RHS[Zeile] += assemblyRHSNeumann(e, A, 2);
                }
            }
        }//close I-Schleife (X-Achse)

        //Kante 10:
        e.resize(2);
        A.resize(2);
        Zeile++;
        if(Dirichlet)
        {
            Matrix.sequential_fill(Zeile, 1.0);
            Matrix.end_of_row();
            RHS[Zeile]= RHSVAL;
        }
        else
        {
            e[0]=Zeile-1 -z; A[0]=5;
            e[1]=Zeile-1; A[1]=1;
            assemblyMatrixRow(e, A, column, value);
            for (int m(0); m<12; ++m)
                Matrix.sequential_fill(column[m], value[m]);
            Matrix.end_of_row();
            RHS[Zeile] = assemblyRHSLoad(e, A, f);
            if(Neumann)
            {
                e[0]=Zeile-1-z; A[0]=2;
                e[1]=Zeile-1; A[1]=1;
                RHS[Zeile] += assemblyRHSNeumann(e, A, 2);
                e[0]=Zeile-z; A[0]=3;
                e[1]=Zeile; A[1]=0;
                RHS[Zeile] += assemblyRHSNeumann(e, A, 3);
            }
        }

        for(int j=1; j<_ny-1; j++)
        {
            //Fläche 5:
            e.resize(4);
            A.resize(4);
            Zeile++;
            if(Dirichlet)
            {
                Matrix.sequential_fill(Zeile, 1.0);
                Matrix.end_of_row();
                RHS[Zeile]= RHSVAL;
            }
            else
            {
                e[0]=Zeile -y-z; A[0]=7;
                e[1]=Zeile -z; A[1]=4;
                e[2]=Zeile; A[2]=0;
                e[3]=Zeile -y; A[3]=3;
                assemblyMatrixRow(e, A, column, value);
                for (int m(0); m<18; ++m)
                    Matrix.sequential_fill(column[m], value[m]);
                Matrix.end_of_row();
                RHS[Zeile] = assemblyRHSLoad(e, A, f);
                if(Neumann)
                {
                    e[0]=Zeile -y-z; A[0]=2;
                    e[1]=Zeile -z; A[1]=3;
                    e[2]=Zeile; A[2]=0;
                    e[3]=Zeile -y; A[3]=1;
                    RHS[Zeile] += assemblyRHSNeumann(e, A, 3);
                }
            }

            //Inneres:
            e.resize(8);
            A.resize(8);
            for(int i=1; i<_nx-1; i++)
            {
                Zeile++;
                //if(Dirichlet)
                //{
                //    Matrix.sequential_fill(Zeile, 1.0);
                //    Matrix.end_of_row();
                //    RHS[Zeile]= RHSVAL;
                //}
                //else
                {
std::cout << Zeile << std::endl;
                    e[0]=Zeile -1-y-z; A[0]=6;
                    e[1]=Zeile -y-z; A[1]=7;
                    e[2]=Zeile -z; A[2]=4;
                    e[3]=Zeile -1-z; A[3]=5;
                    e[4]=Zeile -1-y; A[4]=2;
                    e[5]=Zeile -y; A[5]=3;
                    e[6]=Zeile; A[6]=0;
                    e[7]=Zeile -1; A[7]=1;
                    assemblyMatrixRow(e, A, column, value);
                    //assemblyMatrixRow(std::vector<int>{Zeile-1-y-z, Zeile-y-z, Zeile-z, Zeile-1-z, Zeile-1-y, Zeile-y, Zeile, Zeile-1},
                    //                  std::vector<int>{6, 7, 4, 5, 2, 3, 0, 1}, column, value);
                    for (int m(0); m<27; ++m)
                        Matrix.sequential_fill(column[m], value[m]);
                    Matrix.end_of_row();
                    RHS[Zeile] = assemblyRHSLoad(e, A, f);

                    //Neumann eventuell hinzufuegen
                }
            } //Close I-Schleife (X-Achse)

            //Fläche 6:
            e.resize(4);
            A.resize(4);
            Zeile++;
            if(Dirichlet)
            {
                Matrix.sequential_fill(Zeile, 1.0);
                Matrix.end_of_row();
                RHS[Zeile]= RHSVAL;
            }
            else
            {
                e[0]=Zeile -y-z -1; A[0]=6;
                e[1]=Zeile -z -1; A[1]=5;
                e[2]=Zeile -1; A[2]=1;
                e[3]=Zeile -y -1; A[3]=2;
                assemblyMatrixRow(e, A, column, value);
                for (int m(0); m<18; ++m)
                    Matrix.sequential_fill(column[m], value[m]);
                Matrix.end_of_row();
                RHS[Zeile] = assemblyRHSLoad(e, A, f);
                if(Neumann)
                {
                    e[0]=Zeile -y-z; A[0]=2;
                    e[1]=Zeile -z; A[1]=3;
                    e[2]=Zeile; A[2]=0;
                    e[3]=Zeile -y; A[3]=1;
                    RHS[Zeile] += assemblyRHSNeumann(e, A, 3);
                }
            }
        } //close J-Schleife (Y-Achse)

        //Kante 12:
        e.resize(2);
        A.resize(2);
        Zeile++;
        if(Dirichlet)
        {
            Matrix.sequential_fill(Zeile, 1.0);
            Matrix.end_of_row();
            RHS[Zeile]= RHSVAL;
        }
        else
        {
            e[0]=Zeile-y -z; A[0]=7;
            e[1]=Zeile-y; A[1]=3;
            assemblyMatrixRow(e, A, column, value);
            for (int m(0); m<12; ++m)
                Matrix.sequential_fill(column[m], value[m]);
            Matrix.end_of_row();
            RHS[Zeile] = assemblyRHSLoad(e, A, f);
            if(Neumann)
            {
                e[0]=Zeile-z; A[0]=3;
                e[1]=Zeile; A[1]=0;
                RHS[Zeile] += assemblyRHSNeumann(e, A, 2);
                e[0]=Zeile-z-y; A[0]=2;
                e[1]=Zeile-y; A[1]=1;
                RHS[Zeile] += assemblyRHSNeumann(e, A, 3);
            }
        }

        //Fläche 4
        e.resize(4);
        A.resize(4);
        for(int i=1; i< _nx-1; i++)
        {
            Zeile++;
            if(Dirichlet)
            {
                Matrix.sequential_fill( Zeile, 1.0);
                Matrix.end_of_row();
                RHS[Zeile]= RHSVAL;
            }
            else
            {
                e[0]= Zeile -1-z -y; A[0]=6;
                e[1]= Zeile -z -y; A[1]=7;
                e[2]= Zeile -y; A[2]=3;
                e[3]= Zeile -1 -y; A[3]=2;
                assemblyMatrixRow(e, A, column, value);
                for (int m(0); m<18; ++m)
                    Matrix.sequential_fill(column[m], value[m]);
                Matrix.end_of_row();
                RHS[Zeile] = assemblyRHSLoad(e, A, f);
                if(Neumann)
                {
                    e[0]= Zeile -1-z; A[0]=2;
                    e[1]= Zeile -z; A[1]=3;
                    e[2]= Zeile; A[2]=0;
                    e[3]= Zeile -1; A[3]=1;
                    RHS[Zeile] += assemblyRHSNeumann(e, A, 2);
                }
            }
        }//Close I-Schleife (X-Achse)

        //Kante 11:
        e.resize(2);
        A.resize(2);
        Zeile++;
        if(Dirichlet)
        {
            Matrix.sequential_fill(Zeile, 1.0);
            Matrix.end_of_row();
            RHS[Zeile]= RHSVAL;
        }
        else
        {
            e[0]=Zeile-1-y-z; A[0]=6;
            e[1]=Zeile-1-y; A[1]=2;
            assemblyMatrixRow(e, A, column, value);
            for (int m(0); m<12; ++m)
                Matrix.sequential_fill(column[m], value[m]);
            Matrix.end_of_row();
            RHS[Zeile] = assemblyRHSLoad(e, A, f);
            if(Neumann)
            {
                e[0]=Zeile-1-z; A[0]=2;
                e[1]=Zeile-1; A[1]=1;
                RHS[Zeile] += assemblyRHSNeumann(e, A, 2);
                e[0]=Zeile-z-y; A[0]=2;
                e[1]=Zeile-y; A[1]=1;
                RHS[Zeile] += assemblyRHSNeumann(e, A, 3);
            }
        }
    } //close K-schleife (Z-Achse)

    //Ecke 5
    e.resize(1);
    A.resize(1);
    Zeile++; //Zeile sollte hier (_nz-1)*z sein
    if(Dirichlet)
    {
        Matrix.sequential_fill( Zeile, 1.0);
        Matrix.end_of_row();
        RHS[Zeile]= RHSVAL;
    }
    else
    {
        e[0]= Zeile-z; A[0]=4;
        assemblyMatrixRow(e, A, column, value);
        for (int m(0); m<8; ++m)
            Matrix.sequential_fill(column[m], value[m]);
        Matrix.end_of_row();
        RHS[Zeile] = assemblyRHSLoad(e, A, f);
        if(Neumann)
        {
            e[0]= Zeile; A[0]=0;
            RHS[Zeile] += assemblyRHSNeumann(e, A, 1);
            e[0]= Zeile-z; A[0]=3;
            RHS[Zeile] += assemblyRHSNeumann(e, A, 2);
            e[0]= Zeile-z; A[0]=3;
            RHS[Zeile] += assemblyRHSNeumann(e, A, 3);
        }
    }

    //Kante 4
    e.resize(2);
    A.resize(2);
    for(int i =1; i<_nx-1;i++)
    {
        Zeile++;
        if(Dirichlet)
        {
            Matrix.sequential_fill(Zeile, 1.0);
            Matrix.end_of_row();
            RHS[Zeile]= RHSVAL;
        }
        else
        {
            e[0]=Zeile-z -1; A[0]=5;
            e[1]=Zeile-z; A[1]=4;
            assemblyMatrixRow(e, A, column, value);
            for (int m(0); m<12; ++m)
                Matrix.sequential_fill(column[m], value[m]);
            Matrix.end_of_row();
            RHS[Zeile] = assemblyRHSLoad(e, A, f);
            if(Neumann)
            {
                e[0]=Zeile-1; A[0]=1;
                e[1]=Zeile; A[1]=0;
                RHS[Zeile] += assemblyRHSNeumann(e, A, 1);
                e[0]=Zeile-z -1; A[0]=2;
                e[1]=Zeile-z; A[1]=3;
                RHS[Zeile] += assemblyRHSNeumann(e, A, 2);
            }
        }
    }//Close I-Schleife (X-Achse)

    //Ecke 6:
    e.resize(1);
    A.resize(1);
    Zeile++; //Zeile sollte hier (_nz-1)*z+y-1 sein
    if(Dirichlet)
    {
        Matrix.sequential_fill( Zeile, 1.0);
        Matrix.end_of_row();
        RHS[Zeile]= RHSVAL;
    }
    else
    {
        e[0]= Zeile-1-z; A[0]=5;
        assemblyMatrixRow(e, A, column, value);
        for (int m(0); m<8; ++m)
            Matrix.sequential_fill(column[m], value[m]);
        Matrix.end_of_row();
        RHS[Zeile] = assemblyRHSLoad(e, A, f);
        if(Neumann)
        {
            e[0]= Zeile-1; A[0]=1;
            RHS[Zeile] += assemblyRHSNeumann(e, A, 1);
            e[0]= Zeile-1-z; A[0]=2;
            RHS[Zeile] += assemblyRHSNeumann(e, A, 2);
            e[0]= Zeile-z; A[0]=3;
            RHS[Zeile] += assemblyRHSNeumann(e, A, 3);
        }
    }

    for(int j=1; j< _ny-1; j++)
    {
        //Kante 8
        e.resize(2);
        A.resize(2);
        Zeile++;
        if(Dirichlet)
        {
            Matrix.sequential_fill(Zeile, 1.0);
            Matrix.end_of_row();
            RHS[Zeile]= RHSVAL;
        }
        else
        {
            e[0]=Zeile -z-y; A[0]=7;
            e[1]=Zeile -z; A[1]=4;
            assemblyMatrixRow(e, A, column, value);
            for (int m(0); m<12; ++m)
                Matrix.sequential_fill(column[m], value[m]);
            Matrix.end_of_row();
            RHS[Zeile] = assemblyRHSLoad(e, A, f);
            if(Neumann)
            {
                e[0]=Zeile-y; A[0]=3;
                e[1]=Zeile; A[1]=0;
                RHS[Zeile] += assemblyRHSNeumann(e, A, 1);
                e[0]=Zeile-y-z; A[0]=2;
                e[1]=Zeile-z; A[1]=3;
                RHS[Zeile] += assemblyRHSNeumann(e, A, 3);
            }
        }

        //Fläche 2:
        e.resize(4);
        A.resize(4);
        for(int i=1; i<_nx-1; i++)
        {
            Zeile++;
            if(Dirichlet)
            {
                Matrix.sequential_fill(Zeile, 1.0);
                Matrix.end_of_row();
                RHS[Zeile]= RHSVAL;
            }
            else
            {
                e[0]=Zeile -y-1 -z; A[0]=6;
                e[1]=Zeile -y -z; A[1]=7;
                e[2]=Zeile -z ; A[2]=4;
                e[3]=Zeile -1 -z; A[3]=5;
                assemblyMatrixRow(e, A, column, value);
                for (int m(0); m<18; ++m)
                    Matrix.sequential_fill(column[m], value[m]);
                Matrix.end_of_row();
                RHS[Zeile] = assemblyRHSLoad(e, A, f);
                if(Neumann)
                {
                    e[0]=Zeile -y-1; A[0]=2;
                    e[1]=Zeile -y; A[1]=3;
                    e[2]=Zeile; A[2]=0;
                    e[3]=Zeile -1; A[3]=1;
                    RHS[Zeile] += assemblyRHSNeumann(e, A, 1);
                }
            }
        }//Close I-Schleife (X-Achse)

        //Kante 7
        e.resize(2);
        A.resize(2);
        Zeile++;
        if(Dirichlet)
        {
            Matrix.sequential_fill(Zeile, 1.0);
            Matrix.end_of_row();
            RHS[Zeile]= RHSVAL;
        }
        else
        {
            e[0]=Zeile-1-z -y; A[0]=6;
            e[1]=Zeile-1-z; A[1]=5;
            assemblyMatrixRow(e, A, column, value);
            for (int m(0); m<12; ++m)
                Matrix.sequential_fill(column[m], value[m]);
            Matrix.end_of_row();
            RHS[Zeile] = assemblyRHSLoad(e, A, f);
            if(Neumann)
            {
                e[0]=Zeile-y-1; A[0]=2;
                e[1]=Zeile-1; A[1]=1;
                RHS[Zeile] += assemblyRHSNeumann(e, A, 1);
                e[0]=Zeile-y-z; A[0]=2;
                e[1]=Zeile-z; A[1]=3;
                RHS[Zeile] += assemblyRHSNeumann(e, A, 3);
            }
        }
    }//Close J-Schleife (Y-Achse)

    //Ecke 7
    e.resize(1);
    A.resize(1);
    Zeile++; //Zeile sollte hier (_nx*_ny*_nz)-_nx sein
    if(Dirichlet)
    {
        Matrix.sequential_fill(Zeile, 1.0);
        Matrix.end_of_row();
        RHS[Zeile]= RHSVAL;
    }
    else
    {
        e[0]=Zeile-y-z; A[0]=7;
        assemblyMatrixRow(e, A, column, value);
        for (int m(0); m<8; ++m)
            Matrix.sequential_fill(column[m], value[m]);
        Matrix.end_of_row();
        RHS[Zeile] = assemblyRHSLoad(e, A, f);
        if(Neumann)
        {
            e[0]=Zeile-y; A[0]=3;
            RHS[Zeile] += assemblyRHSNeumann(e, A, 1);
            e[0]=Zeile-z; A[0]=3;
            RHS[Zeile] += assemblyRHSNeumann(e, A, 2);
            e[0]=Zeile-y-z; A[0]=2;
            RHS[Zeile] += assemblyRHSNeumann(e, A, 3);
        }
    }

    //Kante 3
    e.resize(2);
    A.resize(2);
    for(int i=1; i<_nx-1; i++)
    {
        Zeile++;
        if(Dirichlet)
        {
            Matrix.sequential_fill( Zeile, 1.0);
            Matrix.end_of_row();
            RHS[Zeile]= RHSVAL;
        }
        else
        {
            e[0]= Zeile-y-z -1; A[0]=6;
            e[1]= Zeile-y-z; A[1]=7;
            assemblyMatrixRow(e, A, column, value);
            for (int m(0); m<12; ++m)
                Matrix.sequential_fill(column[m], value[m]);
            Matrix.end_of_row();
            RHS[Zeile] = assemblyRHSLoad(e, A, f);
            if(Neumann)
            {
                e[0]= Zeile-1-y; A[0]=2;
                e[1]= Zeile-y; A[1]=3;
                RHS[Zeile] += assemblyRHSNeumann(e, A, 1);
                e[0]= Zeile-1-z; A[0]=2;
                e[1]= Zeile-z; A[1]=3;
                RHS[Zeile] += assemblyRHSNeumann(e, A, 2);
            }
        }
    }//Close I-Schleife (X-Achse)

    //Ecke 8:
    e.resize(1);
    A.resize(1);
    Zeile++; //Zeile sollte hier (_nx*_ny*_nz) sein
    if(Dirichlet)
    {
        Matrix.sequential_fill(Zeile, 1.0);
        Matrix.end_of_row();
        RHS[Zeile]= RHSVAL;
    }
    else
    {
        e[0]=Zeile-z-y-1; A[0]=6;
        assemblyMatrixRow(e, A, column, value);
        for (int m(0); m<8; ++m)
            Matrix.sequential_fill(column[m], value[m]);
        Matrix.end_of_row();
        RHS[Zeile] = assemblyRHSLoad(e, A, f);
        if(Neumann)
        {
            e[0]=Zeile-y-1; A[0]=2;
            RHS[Zeile] += assemblyRHSNeumann(e, A, 1);
            e[0]=Zeile-z-1; A[0]=2;
            RHS[Zeile] += assemblyRHSNeumann(e, A, 2);
            e[0]=Zeile-z-y; A[0]=2;
            RHS[Zeile] += assemblyRHSNeumann(e, A, 3);
        }
    }

    //TODO rhs direkt fuellen (erst wenn alles laeuft)
    for (int i(0); i<_nx*_ny*_nz; ++i)
        rhs.set_global(i, RHS[i]);

}//assemble()

}//namespace Icarus
