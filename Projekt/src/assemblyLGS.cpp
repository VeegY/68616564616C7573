//Ich lasse eine ganze Menge an Redundanz drin, einfach nur für die Lesbarkeit
//Kann man nachher natürlich abändern

int Nx=3;
int Ny=10;
int Nz=5;
int h=1;
int z=Nx*Ny;
int y=Nx;

//TODO: wieder loeschen
#define Dirichlet true
#define Neumann true
//TODO: wieder loeschen


#include "include/assemblyMatrixRow.hpp"
#include "include/assemblyRHSLoad.hpp"
#include "include/assemblyRHSNeumann.hpp"
#include "include/distellpackmatrix.hpp"

namespace Icarus
{

int nomain()
{

    DistEllpackMatrix<double> Matrix(Nx*Ny*Nz);

    int i;
    std::vector<int> e(1);
    std::vector<int> A(1);
    std::vector<double> RHS(Nx*Ny*Nz);

    //Ecke 1
    i=0;
    if(Dirichlet)
    {
        Matrix.prepare_sequential_fill(1);
        Matrix.sequential_fill(i, 1.0);
        Matrix.end_of_row();
        RHS[i]= 999999;
    }
    else
    {
        e[0]=i; A[0]=0;
        Matrix.setZeile(i, assemblyMatrixRow(e, A));
        RHS[i] = assemblyRHSLoad(e, A);
        if(Neumann)
        {
            e[0]=i; A[0]=0;
            RHS[i] += assemblyRHSNeumann(e, A, 1);
            e[0]=i; A[0]=0;
            RHS[i] += assemblyRHSNeumann(e, A, 2);
            e[0]=i; A[0]=0;
            RHS[i] += assemblyRHSNeumann(e, A, 3);
        }
    }

    //Ecke 2
    i=y-1;
    if(Dirichlet)
    {
        Matrix.prepare_sequential_fill(1);
        Matrix.sequential_fill(i, 1.0);
        Matrix.end_of_row();
        RHS[i]= 999999;
    }
    else
    {
        e[0]=i-1; A[0]=1;
        Matrix.setZeile(i, assemblyMatrixRow(e, A));
        RHS[i] = assemblyRHSLoad(e, A);
        if(Neumann)
        {
            e[0]=i-1; A[0]=1;
            RHS[i] += assemblyRHSNeumann(e, A, 1);
            e[0]=i-1; A[0]=1;
            RHS[i] += assemblyRHSNeumann(e, A, 2);
            e[0]=i; A[0]=0;
            RHS[i] += assemblyRHSNeumann(e, A, 3);
        }
    }

    //Ecke 3
    i= z-1;
    if(Dirichlet)
    {
        Matrix.prepare_sequential_fill(1);
        Matrix.sequential_fill(i, 1.0);
        Matrix.end_of_row();
        RHS[i]= 999999;
    }
    else
    {
        e[0]=i-1-y; A[0]=2;
        Matrix.setZeile(i, assemblyMatrixRow(e, A));
        RHS[i] = assemblyRHSLoad(e, A);
        if(Neumann)
        {
            e[0]=i-1-y; A[0]=2;
            RHS[i] += assemblyRHSNeumann(e, A, 1);
            e[0]=i-1; A[0]=1;
            RHS[i] += assemblyRHSNeumann(e, A, 2);
            e[0]=i-y; A[0]=1;
            RHS[i] += assemblyRHSNeumann(e, A, 3);
        }
    }

    //Ecke 4
    i=(Ny-1)*y;
    if(Dirichlet)
    {
        Matrix.prepare_sequential_fill(1);
        Matrix.sequential_fill(i, 1.0);
        Matrix.end_of_row();
        RHS[i]= 999999;
    }
    else
    {
        e[0]=i-y; A[0]=3;
        Matrix.setZeile(i, assemblyMatrixRow(e, A));
        RHS[i] = assemblyRHSLoad(e, A);
        if(Neumann)
        {
            e[0]=i-y; A[0]=3;
            RHS[i] += assemblyRHSNeumann(e, A, 1);
            e[0]=i; A[0]=0;
            RHS[i] += assemblyRHSNeumann(e, A, 2);
            e[0]=i-y; A[0]=1;
            RHS[i] += assemblyRHSNeumann(e, A, 3);
        }
    }

    //Ecke 5
    i=(Nz-1)*z;
    if(Dirichlet)
    {
        Matrix.prepare_sequential_fill(1);
        Matrix.sequential_fill(i, 1.0);
        Matrix.end_of_row();
        RHS[i]= 999999;
    }
    else
    {
        e[0]=i-z; A[0]=4;
        Matrix.setZeile(i, assemblyMatrixRow(e, A));
        RHS[i] = assemblyRHSLoad(e, A);
        if(Neumann)
        {
            e[0]=i; A[0]=0;
            RHS[i] += assemblyRHSNeumann(e, A, 1);
            e[0]=i-z; A[0]=3;
            RHS[i] += assemblyRHSNeumann(e, A, 2);
            e[0]=i-z; A[0]=3;
            RHS[i] += assemblyRHSNeumann(e, A, 3);
        }
    }

    //Ecke 6
    i=(Nz-1)*z+y-1;
    if(Dirichlet)
    {
        Matrix.prepare_sequential_fill(1);
        Matrix.sequential_fill(i, 1.0);
        Matrix.end_of_row();
        RHS[i]= 999999;
    }
    else
    {
        e[0]=i-x-z; A[0]=5;  //TODO TOCHECK: sicher, dass da ein x hinkommt? kompiliert so nicht. vielleicht ist ein y gewollt?
        Matrix.setZeile(i, assemblyMatrixRow(e, A));
        RHS[i] = assemblyRHSLoad(e, A);
        if(Neumann)
        {
            e[0]=i-1; A[0]=1;
            RHS[i] += assemblyRHSNeumann(e, A, 1);
            e[0]=i-1-z; A[0]=2;
            RHS[i] += assemblyRHSNeumann(e, A, 2);
            e[0]=i-z; A[0]=3;
            RHS[i] += assemblyRHSNeumann(e, A, 3);
        }
    }

    //Ecke 7
    i=Nx*Ny*Nz -1;
    if(Dirichlet)
    {
        Matrix.prepare_sequential_fill(1);
        Matrix.sequential_fill(i, 1.0);
        Matrix.end_of_row();
        RHS[i]= 999999;
    }
    else
    {
        e[0]=i-1-y-z; A[0]=6;
        Matrix.setZeile(i, assemblyMatrixRow(e, A));
        RHS[i] = assemblyRHSLoad(e, A);
        if(Neumann)
        {
            e[0]=i-1-y; A[0]=2;
            RHS[i] += assemblyRHSNeumann(e, A, 1);
            e[0]=i-1-z; A[0]=2;
            RHS[i] += assemblyRHSNeumann(e, A, 2);
            e[0]=i-y-z; A[0]=2;
            RHS[i] += assemblyRHSNeumann(e, A, 3);
        }
    }

    //Ecke 8
    i=(z-1)*y;
    if(Dirichlet)
    {
        Matrix.prepare_sequential_fill(1);
        Matrix.sequential_fill(i, 1.0);
        Matrix.end_of_row();
        RHS[i]= 999999;
    }
    else
    {
        e[0]=i-y-z; A[0]=7;
        Matrix.setZeile(i, assemblyMatrixRow(e, A));
        RHS[i] = assemblyRHSLoad(e, A);
        if(Neumann)
        {
            e[0]=i-y; A[0]=3;
            RHS[i] += assemblyRHSNeumann(e, A, 1);
            e[0]=i-z; A[0]=3;
            RHS[i] += assemblyRHSNeumann(e, A, 2);
            e[0]=i-z-y; A[0]=2;
            RHS[i] += assemblyRHSNeumann(e, A, 3);
        }
    }

    e.resize(2);
    A.resize(2);

    //Kante 1:
    for(int j=1; j<Nx-1;j++)
    {
        i=j;
        if(Dirichlet)
        {
            Matrix.prepare_sequential_fill(1);
            Matrix.sequential_fill(i, 1.0);
            Matrix.end_of_row();
            RHS[i]= 999999;
        }
        else
        {
            e[0]=i-1; A[0]=1;
            e[1]=i; A[1]=0;
            Matrix.setZeile(i, assemblyMatrixRow(e, A));
            RHS[i] = assemblyRHSLoad(e, A);
            if(Neumann)
            {
                e[0]=i-1; A[0]=1;
                e[1]=i; A[1]=0;
                RHS[i] += assemblyRHSNeumann(e, A, 1);
                e[0]=i-1; A[0]=1;
                e[1]=i; A[1]=0;
                RHS[i] += assemblyRHSNeumann(e, A, 2);
            }
        }
    }

    //Kante 2:
    for(int j=1; j<Nx-1; j++)
    {
        i=(Ny-1)*y + j;
        if(Dirichlet)
        {
            Matrix.prepare_sequential_fill(1);
            Matrix.sequential_fill(i, 1.0);
            Matrix.end_of_row();
            RHS[i]= 999999;
        }
        else
        {
            e[0]=i-y -1; A[0]=2;
            e[1]=i-y; A[1]=3;
            Matrix.setZeile(i, assemblyMatrixRow(e, A));
            RHS[i] = assemblyRHSLoad(e, A);
            if(Neumann)
            {
                e[0]=i-y-1; A[0]=2;
                e[1]=i-y; A[1]=3;
                RHS[i] += assemblyRHSNeumann(e, A, 1);
                e[0]=i-1; A[0]=1;
                e[1]=i; A[1]=0;
                RHS[i] += assemblyRHSNeumann(e, A, 2);
            }
        }
    }

    //Kante 3:
    for(int j=1; j<Nx-1; j++)
    {
        i=(z-1)*y + j;
        if(Dirichlet)
        {
            Matrix.prepare_sequential_fill(1);
            Matrix.sequential_fill(i, 1.0);
            Matrix.end_of_row();
            RHS[i]= 999999;
        }
        else
        {
            e[0]=i-y-z -1; A[0]=6;
            e[1]=i-y-z; A[1]=7;
            Matrix.setZeile(i, assemblyMatrixRow(e, A));
            RHS[i] = assemblyRHSLoad(e, A);
            if(Neumann)
            {
                e[0]=i-1-y; A[0]=2;
                e[1]=i-y; A[1]=3;
                RHS[i] += assemblyRHSNeumann(e, A, 1);
                e[0]=i-1-z; A[0]=2;
                e[1]=i-z; A[1]=3;
                RHS[i] += assemblyRHSNeumann(e, A, 2);
            }
        }
    }

    //Kante 4:
    for(int j=1; j<Nx-1; j++)
    {
       i=(Nz-1)*z + j;
        if(Dirichlet)
        {
            Matrix.prepare_sequential_fill(1);
            Matrix.sequential_fill(i, 1.0);
            Matrix.end_of_row();
            RHS[i]= 999999;
        }
        else
        {
            e[0]=i-z -1; A[0]=5;
            e[1]=i-z; A[1]=4;
            Matrix.setZeile(i, assemblyMatrixRow(e, A));
            RHS[i] = assemblyRHSLoad(e, A);
            if(Neumann)
            {
                e[0]=i-1; A[0]=1;
                e[1]=i; A[1]=0;
                RHS[i] += assemblyRHSNeumann(e, A, 1);
                e[0]=i-z -1; A[0]=2;
                e[1]=i-z; A[1]=3;
                RHS[i] += assemblyRHSNeumann(e, A, 2);
            }
        }
    }

    //Kante 5:
    for(int j=1; j<Ny-1; j++)
    {
        i= j*y;
        if(Dirichlet)
        {
            Matrix.prepare_sequential_fill(1);
            Matrix.sequential_fill(i, 1.0);
            Matrix.end_of_row();
            RHS[i]= 999999;
        }
        else
        {
            e[0]=i-y; A[0]=3;
            e[1]=i; A[1]=0;
            Matrix.setZeile(i, assemblyMatrixRow(e, A));
            RHS[i] = assemblyRHSLoad(e, A);
            if(Neumann)
            {
                e[0]=i-y; A[0]=3;
                e[1]=i; A[1]=0;
                RHS[i] += assemblyRHSNeumann(e, A, 1);
                e[0]=i-y; A[0]=1;
                e[1]=i; A[1]=0;
                RHS[i] += assemblyRHSNeumann(e, A, 3);
            }
        }
    }

    //Kante 6:
    for(int j=1; j<Ny-1; j++)
    {
        i=(y-1) + j*y;
        if(Dirichlet)
        {
            Matrix.prepare_sequential_fill(1);
            Matrix.sequential_fill(i, 1.0);
            Matrix.end_of_row();
            RHS[i]= 999999;
        }
        else
        {
            e[0]=i-1 -y; A[0]=2;
            e[1]=i-1; A[1]=1;
            Matrix.setZeile(i, assemblyMatrixRow(e, A));
            RHS[i] = assemblyRHSLoad(e, A);
            if(Neumann)
            {
                e[0]=i-y-1; A[0]=2;
                e[1]=i-1; A[1]=1;
                RHS[i] += assemblyRHSNeumann(e, A, 1);
                e[0]=i-y; A[0]=1;
                e[1]=i; A[1]=0;
                RHS[i] += assemblyRHSNeumann(e, A, 3);
            }
        }
    }

    //Kante 7:
    for(int j=1; j<Ny-1; j++)
    {
        i=(Nz-1)*z+y-1 + j*y;
        if(Dirichlet)
        {
            Matrix.prepare_sequential_fill(1);
            Matrix.sequential_fill(i, 1.0);
            Matrix.end_of_row();
            RHS[i]= 999999;
        }
        else
        {
            e[0]=i-1-z -y; A[0]=6;
            e[1]=i-1-z; A[1]=5;
            Matrix.setZeile(i, assemblyMatrixRow(e, A));
            RHS[i] = assemblyRHSLoad(e, A);
            if(Neumann)
            {
                e[0]=i-y-1; A[0]=2;
                e[1]=i-1; A[1]=1;
                RHS[i] += assemblyRHSNeumann(e, A, 1);
                e[0]=i-y-z; A[0]=2;
                e[1]=i-z; A[1]=3;
                RHS[i] += assemblyRHSNeumann(e, A, 3);
            }
        }
    }

    //Kante 8:
    for(int j=1; j<Nz-1; j++)
    {
        i= j*z;
        if(Dirichlet)
        {
            Matrix.prepare_sequential_fill(1);
            Matrix.sequential_fill(i, 1.0);
            Matrix.end_of_row();
            RHS[i]= 999999;
        }
        else
        {
            e[0]=i -z; A[0]=4;
            e[1]=i; A[1]=0;
            Matrix.setZeile(i, assemblyMatrixRow(e, A));
            RHS[i] = assemblyRHSLoad(e, A);
            if(Neumann)
            {
                e[0]=i-y; A[0]=3;
                e[1]=i; A[1]=0;
                RHS[i] += assemblyRHSNeumann(e, A, 1);
                e[0]=i-y-z; A[0]=2;
                e[1]=i-z; A[1]=3;
                RHS[i] += assemblyRHSNeumann(e, A, 3);
            }
        }
    }

    //Kante 9:
    for(int j=1; j<Ny-1; j++)
    {
        i=(Nz-1)*z + j*y;
        if(Dirichlet)
        {
            Matrix.prepare_sequential_fill(1);
            Matrix.sequential_fill(i, 1.0);
            Matrix.end_of_row();
            RHS[i]= 999999;
        }
        else
        {
            e[0]=i-z -y; A[0]=7;
            e[1]=i-z; A[1]=4;
            Matrix.setZeile(i, assemblyMatrixRow(e, A));
            RHS[i] = assemblyRHSLoad(e, A);
            if(Neumann)
            {
                e[0]=i-z; A[0]=3;
                e[1]=i; A[1]=0;
                RHS[i] += assemblyRHSNeumann(e, A, 2);
                e[0]=i-z; A[0]=3;
                e[1]=i; A[1]=0;
                RHS[i] += assemblyRHSNeumann(e, A, 3);
            }
        }
    }

    //Kante 10:
    for(int j=1; j<Ny-1; j++)
    {
        i=y-1 + j*z;
        if(Dirichlet)
        {
            Matrix.prepare_sequential_fill(1);
            Matrix.sequential_fill(i, 1.0);
            Matrix.end_of_row();
            RHS[i]= 999999;
        }
        else
        {
            e[0]=i-1 -z; A[0]=5;
            e[1]=i-1; A[1]=1;
            Matrix.setZeile(i, assemblyMatrixRow(e, A));
            RHS[i] = assemblyRHSLoad(e, A);
            if(Neumann)
            {
                e[0]=i-1-z; A[0]=2;
                e[1]=i-1; A[1]=1;
                RHS[i] += assemblyRHSNeumann(e, A, 2);
                e[0]=i-z; A[0]=3;
                e[1]=i; A[1]=0;
                RHS[i] += assemblyRHSNeumann(e, A, 3);
            }
        }
    }

    //Kante 11:
    for(int j=1; j<Ny-1; j++)
    {
        i=y-1 + j*z;
        if(Dirichlet)
        {
            Matrix.prepare_sequential_fill(1);
            Matrix.sequential_fill(i, 1.0);
            Matrix.end_of_row();
            RHS[i]= 999999;
        }
        else
        {
            e[0]=i-1 -z; A[0]=5;
            e[1]=i-1; A[1]=1;
            Matrix.setZeile(i, assemblyMatrixRow(e, A));
            RHS[i] = assemblyRHSLoad(e, A);
            if(Neumann)
            {
                e[0]=i-1-z; A[0]=2;
                e[1]=i-1; A[1]=1;
                RHS[i] += assemblyRHSNeumann(e, A, 2);
                e[0]=i-z-y; A[0]=2;
                e[1]=i-y; A[1]=1;
                RHS[i] += assemblyRHSNeumann(e, A, 3);
            }
        }
    }

    //Kante 12:
    for(int j=1; j<Ny-1; j++)
    {
        i=z-y + j*z;
        if(Dirichlet)
        {
            Matrix.prepare_sequential_fill(1);
            Matrix.sequential_fill(i, 1.0);
            Matrix.end_of_row();
            RHS[i]= 999999;
        }
        else
        {
            e[0]=i-y -z; A[0]=7;
            e[1]=i-y; A[1]=3;
            Matrix.setZeile(i, assemblyMatrixRow(e, A));
            RHS[i] = assemblyRHSLoad(e, A);
            if(Neumann)
            {
                e[0]=i-z; A[0]=3;
                e[1]=i; A[1]=0;
                RHS[i] += assemblyRHSNeumann(e, A, 2);
                e[0]=i-z-y; A[0]=2;
                e[1]=i-y; A[1]=1;
                RHS[i] += assemblyRHSNeumann(e, A, 3);
            }
        }
    }

    e.resize(4);
    A.resize(4);

    //Fläche 1:
    for(int j=1; j<Nx-1; j++)
    {
        for(int k=1; k<Ny-1; k++)
        {
            i = j + k*y;
            if(Dirichlet)
            {
                Matrix.prepare_sequential_fill(1);
                Matrix.sequential_fill(i, 1.0);
                Matrix.end_of_row();
                RHS[i]= 999999;
            }
            else
            {
                e[0]=i -y-1; A[0]=2;
                e[1]=i -y; A[1]=3;
                e[2]=i; A[2]=0;
                e[3]=i -1; A[3]=1;
                Matrix.setZeile(i, assemblyMatrixRow(e, A));
                RHS[i] = assemblyRHSLoad(e, A);
                if(Neumann)
                {
                    e[0]=i -y-1; A[0]=2;
                    e[1]=i -y; A[1]=3;
                    e[2]=i; A[2]=0;
                    e[3]=i -1; A[3]=1;
                    RHS[i] += assemblyRHSNeumann(e, A, 1);
                }
            }
       }
    }

    //Fläche 2:
    for(int j=1; j<Nx-1; j++)
    {
        for(int k=1; k<Ny-1; k++)
        {
            i = (Nz-1)*z + j + k*y;
            if(Dirichlet)
            {
                Matrix.prepare_sequential_fill(1);
                Matrix.sequential_fill(i, 1.0);
                Matrix.end_of_row();
                RHS[i]= 999999;
            }
            else
            {
                e[0]=i -y-1 -z; A[0]=6;
                e[1]=i -y -z; A[1]=7;
                e[2]=i -z ; A[2]=4;
                e[3]=i -1 -z; A[3]=5;
                Matrix.setZeile(i, assemblyMatrixRow(e, A));
                RHS[i] = assemblyRHSLoad(e, A);
                if(Neumann)
                {
                    e[0]=i -y-1; A[0]=2;
                    e[1]=i -y; A[1]=3;
                    e[2]=i; A[2]=0;
                    e[3]=i -1; A[3]=1;
                    RHS[i] += assemblyRHSNeumann(e, A, 1);
                }
            }
       }
    }

    //Fläche 3:
    for(int j=1; j<Nx-1; j++)
    {
        for(int k=1; k<Nz-1; k++)
        {
            i = j + k*z;
            if(Dirichlet)
            {
                Matrix.prepare_sequential_fill(1);
                Matrix.sequential_fill(i, 1.0);
                Matrix.end_of_row();
                RHS[i]= 999999;
            }
            else
            {
                e[0]=i -1-z; A[0]=5;
                e[1]=i -z; A[1]=4;
                e[2]=i; A[2]=0;
                e[3]=i -1; A[3]=1;
                Matrix.setZeile(i, assemblyMatrixRow(e, A));
                RHS[i] = assemblyRHSLoad(e, A);
                if(Neumann)
                {
                    e[0]=i -1-z; A[0]=2;
                    e[1]=i -z; A[1]=3;
                    e[2]=i; A[2]=0;
                    e[3]=i -1; A[3]=1;
                    RHS[i] += assemblyRHSNeumann(e, A, 2);
                }
            }
       }
    }

    //Fläche 4:
    for(int j=1; j<Nx-1; j++)
    {
        for(int k=1; k<Nz-1; k++)
        {
            i = (Ny-1)*y + j + k*z;
            if(Dirichlet)
            {
                Matrix.prepare_sequential_fill(1);
                Matrix.sequential_fill(i, 1.0);
                Matrix.end_of_row();
                RHS[i]= 999999;
            }
            else
            {
                e[0]=i -1-z -y; A[0]=6;
                e[1]=i -z -y; A[1]=7;
                e[2]=i -y; A[2]=3;
                e[3]=i -1 -y; A[3]=2;
                Matrix.setZeile(i, assemblyMatrixRow(e, A));
                RHS[i] = assemblyRHSLoad(e, A);
                if(Neumann)
                {
                    e[0]=i -1-z; A[0]=2;
                    e[1]=i -z; A[1]=3;
                    e[2]=i; A[2]=0;
                    e[3]=i -1; A[3]=1;
                    RHS[i] += assemblyRHSNeumann(e, A, 2);
                }
            }
       }
    }

    //Fläche 5:
    for(int j=1; j<Ny-1; j++)
    {
        for(int k=1; k<Nz-1; k++)
        {
            i = j*y + k*z;
            if(Dirichlet)
            {
                Matrix.prepare_sequential_fill(1);
                Matrix.sequential_fill(i, 1.0);
                Matrix.end_of_row();
                RHS[i]= 999999;
            }
            else
            {
                e[0]=i -y-z; A[0]=7;
                e[1]=i -z; A[1]=4;
                e[2]=i; A[2]=0;
                e[3]=i -y; A[3]=3;
                Matrix.setZeile(i, assemblyMatrixRow(e, A));
                RHS[i] = assemblyRHSLoad(e, A);
                if(Neumann)
                {
                    e[0]=i -y-z; A[0]=2;
                    e[1]=i -z; A[1]=3;
                    e[2]=i; A[2]=0;
                    e[3]=i -y; A[3]=1;
                    RHS[i] += assemblyRHSNeumann(e, A, 3);
                }
            }
       }
    }

    //Fläche 6:
    for(int j=1; j<Ny-1; j++)
    {
        for(int k=1; k<Nz-1; k++)
        {
            i = (Nx-1) + j*y + k*z;
            if(Dirichlet)
            {
                Matrix.prepare_sequential_fill(1);
                Matrix.sequential_fill(i, 1.0);
                Matrix.end_of_row();
                RHS[i]= 999999;
            }
            else
            {
                e[0]=i -y-z -1; A[0]=6;
                e[1]=i -z -1; A[1]=5;
                e[2]=i -1; A[2]=1;
                e[3]=i -y -1; A[3]=2;
                Matrix.setZeile(i, assemblyMatrixRow(e, A));
                RHS[i] = assemblyRHSLoad(e, A);
                if(Neumann)
                {
                    e[0]=i -y-z; A[0]=2;
                    e[1]=i -z; A[1]=3;
                    e[2]=i; A[2]=0;
                    e[3]=i -y; A[3]=1;
                    RHS[i] += assemblyRHSNeumann(e, A, 3);
                }
            }
       }
    }

    e.resize(8);
    A.resize(8);

    //Inneres 1:
    for(int j=1; j<Nx-1; j++)
    {
        for(int k=1; k<Ny-1; k++)
        {
            for(int l=1; l<Nz-1; l++)
            {
                i= j + k*y + l*z;
                if(Dirichlet)
                {
                    Matrix.prepare_sequential_fill(1);
                    Matrix.sequential_fill(i, 1.0);
                    Matrix.end_of_row();
                    RHS[i]= 999999;
                }
                else
                {
                    e[0]=i -1-y-z; A[0]=6;
                    e[1]=i -y-z; A[1]=7;
                    e[2]=i -z; A[2]=4;
                    e[3]=i -1-z; A[3]=5;
                    e[4]=i -1-y; A[4]=2;
                    e[5]=i -y; A[5]=3;
                    e[6]=i; A[6]=0;
                    e[7]=i -1; A[7]=1;
                    Matrix.setZeile(i, assemblyMatrixRow(e, A));
                    RHS[i] = assemblyRHSLoad(e, A);

                    //Neumann eventuell hinzufuegen
                }
            }
        }
    }

    return 0;

}//nomain()

}//namespace Icarus
