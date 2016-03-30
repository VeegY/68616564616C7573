#include <iostream>


#include "assemble.hpp"

namespace Icarus
{
	template<typename Scalar>
	void assemble_row(
		DistEllpackMatrix<Scalar>& A,
		SlicedVector<Scalar>& rhs,
		typename ScalarTraits<Scalar>::RealType h,
		std::vector<char>& types, int Nx, int Ny, size_t vtx_global, Scalar rhs_val)
	{
		const size_t fron = A.first_row_on_node();
		const size_t vtx_local = vtx_global - fron;

		// Typ boundary or obvtx_globalect
		if (types[vtx_global] == 'b' || types[vtx_global] == 'o')
		{
			A.sequential_fill(vtx_global, 1.0);
			// homog. dirichlet
			rhs.set_local(vtx_local, 0.0);
		}
		// Typ freier knoten
		else
		{
			// nachbarn (x+,x-,y+,y-,z+,z-)
			const size_t nn[6] = {
				vtx_global + 1, vtx_global - 1,
				vtx_global + Nx, vtx_global - Nx,
				vtx_global + Nx*Ny, vtx_global - Nx*Ny };

			A.sequential_fill(vtx_global, -6.0);
			for (int i = 0; i < 6; i++) A.sequential_fill(nn[i], 1.0);

			// rechte seite
			rhs.set_local(vtx_local, rhs_val);
		}
		A.end_of_row();
	}


	template<typename Scalar>
	std::pair < DistEllpackMatrix<Scalar>,
		SlicedVector < Scalar >>
		assemble(std::vector<char>& disc_points,
		typename ScalarTraits<Scalar>::RealType h,
		int Nx, int Ny, int Nz,
		std::function<Scalar(size_t, size_t, size_t)> rhs_func)
	{
			const size_t N = Nx*Ny*Nz;
			DistEllpackMatrix<Scalar> A(N);
			SlicedVector<Scalar> rhs(N);

			const size_t fron = A.first_row_on_node();
			size_t fron_x, fron_y, fron_z;
			deflatten_3d(fron, Nx, Ny, fron_x, fron_y, fron_z);

			const size_t lron = fron + A.get_dim_local() - 1;
			size_t lron_x, lron_y, lron_z;
			deflatten_3d(lron, Nx, Ny, lron_x, lron_y, lron_z);

			// 7 punkte stern
			const unsigned max_nnz_per_line = 7;
			A.prepare_sequential_fill(max_nnz_per_line);
			for (size_t z = fron_z; z <= lron_z; z++)
			{
				size_t ymin = (z == fron_z) ? fron_y : 0;
				size_t ymax = (z == lron_z) ? lron_y : Ny - 1;
				for (size_t y = ymin; y <= ymax; y++)
				{
					size_t xmin = (z == fron_z && y == fron_y) ? fron_x : 0;
					size_t xmax = (z == lron_z && y == lron_y) ? lron_x : Nx - 1;
					for (size_t x = xmin; x <= xmax; x++)
					{
						const size_t index = x + y*Nx + z*Nx*Ny;
						assemble_row(A, rhs, h, disc_points, Nx, Ny, index, rhs_func(x, y, z));
					}
				}
			}
			return{ A, rhs };
		}


	template<typename Scalar>
	std::pair < DistEllpackMatrix<Scalar>,
		SlicedVector < Scalar >>
		assemble_neumann(size_t nx, size_t ny, size_t nz,
		typename ScalarTraits<Scalar>::RealType h,
		std::function<Scalar(size_t)> bdry)
	{
			const size_t N = nx*ny*nz;
			DistEllpackMatrix<Scalar> A(N);
			SlicedVector<Scalar> rhs(N);

			size_t fron = A.first_row_on_node();
			size_t lron = fron + A.get_dim_local() - 1;

			A.prepare_sequential_fill(7);

			for (size_t vtx_global = fron; vtx_global <= lron; vtx_global++)
			{
				std::vector<int> index = { 0, 0, 0, 0, 0, 0, 0 };
				std::vector<Scalar> wert = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

				//vereinfacht den Zugriff auf die Array-Elemente
				int i = vtx_global + 1;

				//Überprüfung der Position
				if (i <= nx*ny) //Boden
				{
					if (i == 1) //vorderer unterer linker Eckpunkt
					{
						//Dirichlet Randwert in der vorderen unteren linken Ecke
						A.sequential_fill(0, 1);

						index[0] = vtx_global;
						index[1] = vtx_global + 1;
						index[2] = vtx_global + 2;
						index[3] = vtx_global + nx;
						index[4] = vtx_global + 2 * nx;
						index[5] = vtx_global + nx*ny;
						index[6] = vtx_global + 2 * nx*ny;

						//zentraler Differenzenquotient in keine Richtung möglich

						wert[0] = 3.0 * 11.0 / 38.0;
						wert[1] = -28.0 / 38.0;
						wert[2] = 17.0 / 38.0;
						wert[3] = -28.0 / 38.0;
						wert[4] = 17.0 / 38.0;
						wert[5] = -28.0 / 38.0;
						wert[6] = 17.0 / 38.0;

						//zeile[vtx_global] = 3.0*11.0/38.0;
						//zeile[vtx_global+1] = -28.0/38.0;
						//zeile[vtx_global+2] = 17.0/38.0;
						//zeile[vtx_global+nx] = -28.0/38.0;
						//zeile[vtx_global+2*nx] = 17.0/38.0;
						//zeile[vtx_global+nx*ny] = -28.0/38.0;
						//zeile[vtx_global+2*nx*ny] = 17.0/38.0;

						//NeumannRB, Normalenvektor ist (1/sqrt(3),1/sqrt(3),1/sqrt(3))
						//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten

						wert[0] += 3.0 * 1.0 / sqrt(3.0) * 3.0 / 2.0 * h;
						wert[1] += 1.0 / sqrt(3.0)*(-h) / 2.0;
						wert[2] += 1.0 / sqrt(3.0) * 2.0 * h;
						wert[3] += 1.0 / sqrt(3.0)*(-h) / 2.0;
						wert[4] += 1.0 / sqrt(3.0) * 2.0 * h;
						wert[5] += 1.0 / sqrt(3.0)*(-h) / 2.0;
						wert[6] += 1.0 / sqrt(3.0) * 2.0 * h;

						//zeile[vtx_global] += 3.0*1.0/sqrt(3.0)*3.0/2.0*h;
						//zeile[vtx_global+1] += 1.0/sqrt(3.0)*(-h)/2.0;
						//zeile[vtx_global+2] += 1.0/sqrt(3.0) *2.0*h;
						//zeile[vtx_global+nx] += 1.0/sqrt(3.0)*(-h)/2.0;
						//zeile[vtx_global+2*nx] += 1.0/sqrt(3.0) *2.0*h;
						//zeile[vtx_global+nx*ny] += 1.0/sqrt(3.0)*(-h)/2.0;
						//zeile[vtx_global+2*nx*ny] += 1.0/sqrt(3.0) *2.0*h;

					}

					else if (i == nx) //vorderer unterer rechter Eckpunkt
					{
						index[0] = vtx_global;
						index[1] = vtx_global - 1;
						index[2] = vtx_global - 2;
						index[3] = vtx_global + nx;
						index[4] = vtx_global + 2 * nx;
						index[5] = vtx_global + nx*ny;
						index[6] = vtx_global + 2 * nx*ny;

						//zentraler Differenzenquotient in keine Richtung möglich

						wert[0] = 3.0 * 11.0 / 38.0;
						wert[1] = -28.0 / 38.0;
						wert[2] = 17.0 / 38.0;
						wert[3] = -28.0 / 38.0;
						wert[4] = 17.0 / 38.0;
						wert[5] = -28.0 / 38.0;
						wert[6] = 17.0 / 38.0;


						//zeile[vtx_global] = 3.0*11.0/38.0;
						//zeile[vtx_global-1] = -28.0/38.0;
						//zeile[vtx_global-2] = 17.0/38.0;
						//zeile[vtx_global+nx] = -28.0/38.0;
						//zeile[vtx_global+2*nx] = 17.0/38.0;
						//zeile[vtx_global+nx*ny] = -28.0/38.0;
						//zeile[vtx_global+2*nx*ny] = 17.0/38.0;

						//NeumannRB, Normalenvektor ist (-1/sqrt(3),1/sqrt(3),1/sqrt(3))
						//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten

						wert[0] += 1.0 * 1.0 / sqrt(3.0) * 3.0 / 2.0 * h;
						wert[1] += (-1.0) * 1.0 / sqrt(3.0)*(-h) / 2.0;
						wert[2] += (-1.0) * 1.0 / sqrt(3.0) * 2.0 * h;
						wert[3] += 1.0 / sqrt(3.0)*(-h) / 2.0;
						wert[4] += 1.0 / sqrt(3.0) * 2.0 * h;
						wert[5] += 1.0 / sqrt(3.0)*(-h) / 2.0;
						wert[6] += 1.0 / sqrt(3.0) * 2.0 * h;

						//zeile[vtx_global] += 1.0*1.0/sqrt(3.0)*3.0/2.0*h;
						//zeile[vtx_global-1] += (-1.0)*1.0/sqrt(3.0)*(-h)/2.0;
						//zeile[vtx_global-2] += (-1.0)*1.0/sqrt(3.0) *2.0*h;
						//zeile[vtx_global+nx] += 1.0/sqrt(3.0)*(-h)/2.0;
						//zeile[vtx_global+2*nx] += 1.0/sqrt(3.0) *2.0*h;
						//zeile[vtx_global+nx*ny] += 1.0/sqrt(3.0)*(-h)/2.0;
						//zeile[vtx_global+2*nx*ny] += 1.0/sqrt(3.0) *2.0*h;
					}

					else if (i == nx*(ny - 1) + 1) //hinterer unterer linker Eckpunkt
					{
						index[0] = vtx_global;
						index[1] = vtx_global + 1;
						index[2] = vtx_global + 2;
						index[3] = vtx_global - nx;
						index[4] = vtx_global - 2 * nx;
						index[5] = vtx_global + nx*ny;
						index[6] = vtx_global + 2 * nx*ny;

						//zentraler Differenzenquotient in keine Richtung möglich

						wert[0] = 3.0 * 11.0 / 38.0;
						wert[1] = -28.0 / 38.0;
						wert[2] = 17.0 / 38.0;
						wert[3] = -28.0 / 38.0;
						wert[4] = 17.0 / 38.0;
						wert[5] = -28.0 / 38.0;
						wert[6] = 17.0 / 38.0;

						//zeile[vtx_global] = 3.0*11.0/38.0;
						//zeile[vtx_global+1] = -28.0/38.0;
						//zeile[vtx_global+2] = 17.0/38.0;
						//zeile[vtx_global-nx] = -28.0/38.0;
						//zeile[vtx_global-2*nx] = 17.0/38.0;
						//zeile[vtx_global+nx*ny] = -28.0/38.0;
						//zeile[vtx_global+2*nx*ny] = 17.0/38.0;

						//NeumannRB, Normalenvektor ist (1/sqrt(3),-1/sqrt(3),1/sqrt(3))
						//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten

						wert[0] += 1.0 * 1.0 / sqrt(3.0) * 3.0 / 2.0 * h;
						wert[1] += 1.0 / sqrt(3.0)*(-h) / 2.0;
						wert[2] += 1.0 / sqrt(3.0) * 2.0 * h;
						wert[3] += (-1.0) * 1.0 / sqrt(3.0)*(-h) / 2.0;
						wert[4] += (-1.0) * 1.0 / sqrt(3.0) * 2.0 * h;
						wert[5] += 1.0 / sqrt(3.0)*(-h) / 2.0;
						wert[6] += 1.0 / sqrt(3.0) * 2.0 * h;

						//zeile[vtx_global] += 1.0*1.0/sqrt(3.0)*3.0/2.0*h;
						//zeile[vtx_global+1] += 1.0/sqrt(3.0)*(-h)/2.0;
						//zeile[vtx_global+2] += 1.0/sqrt(3.0) *2.0*h;
						//zeile[vtx_global-nx] += (-1.0)*1.0/sqrt(3.0)*(-h)/2.0;
						//zeile[vtx_global-2*nx] += (-1.0)*1.0/sqrt(3.0) *2.0*h;
						//zeile[vtx_global+nx*ny] += 1.0/sqrt(3.0)*(-h)/2.0;
						//zeile[vtx_global+2*nx*ny] += 1.0/sqrt(3.0) *2.0*h;
					}

					else if (i == nx*ny) //hinterer  unterer rechter Eckpunkt
					{
						index[0] = vtx_global;
						index[1] = vtx_global - 1;
						index[2] = vtx_global - 2;
						index[3] = vtx_global - nx;
						index[4] = vtx_global - 2 * nx;
						index[5] = vtx_global + nx*ny;
						index[6] = vtx_global + 2 * nx*ny;

						//zentraler Differenzenquotient in keine Richtung möglich

						wert[0] = 3.0 * 11.0 / 38.0;
						wert[1] = -28.0 / 38.0;
						wert[2] = 17.0 / 38.0;
						wert[3] = -28.0 / 38.0;
						wert[4] = 17.0 / 38.0;
						wert[5] = -28.0 / 38.0;
						wert[6] = 17.0 / 38.0;

						//zeile[vtx_global] = 3.0*11.0/38.0;
						//zeile[vtx_global-1] = -28.0/38.0;
						//zeile[vtx_global-2] = 17.0/38.0;
						//zeile[vtx_global-nx] = -28.0/38.0;
						//zeile[vtx_global-2*nx] = 17.0/38.0;
						//zeile[vtx_global+nx*ny] = -28.0/38.0;
						//zeile[vtx_global+2*nx*ny] = 17.0/38.0;

						//NeumannRB, Normalenvektor ist (-1/sqrt(3),-1/sqrt(3),1/sqrt(3))
						//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten

						wert[0] += (-1.0) * 1.0 / sqrt(3.0) * 3.0 / 2.0 * h;
						wert[1] += (-1.0) * 1.0 / sqrt(3.0)*(-h) / 2.0;
						wert[2] += (-1.0) * 1.0 / sqrt(3.0) * 2.0 * h;
						wert[3] += (-1.0) * 1.0 / sqrt(3.0)*(-h) / 2.0;
						wert[4] += (-1.0) * 1.0 / sqrt(3.0) * 2.0 * h;
						wert[5] += 1.0 / sqrt(3.0)*(-h) / 2.0;
						wert[6] += 1.0 / sqrt(3.0) * 2.0 * h;

						//zeile[vtx_global] += (-1.0)*1.0/sqrt(3.0)*3.0/2.0*h;
						//zeile[vtx_global-1] += (-1.0)*1.0/sqrt(3.0)*(-h)/2.0;
						//zeile[vtx_global-2] += (-1.0)*1.0/sqrt(3.0) *2.0*h;
						//zeile[vtx_global-nx] += (-1.0)*1.0/sqrt(3.0)*(-h)/2.0;
						//zeile[vtx_global-2*nx] += (-1.0)*1.0/sqrt(3.0) *2.0*h;
						//zeile[vtx_global+nx*ny] += 1.0/sqrt(3.0)*(-h)/2.0;
						//zeile[vtx_global+2*nx*ny] += 1.0/sqrt(3.0) *2.0*h;
					}

					else if (i < nx) //vordere untere Kante ohne Eckpunkt
					{
						index[0] = vtx_global;
						index[1] = vtx_global + 1;
						index[2] = vtx_global - 1;
						index[3] = vtx_global + nx;
						index[4] = vtx_global + 2 * nx;
						index[5] = vtx_global + nx*ny;
						index[6] = vtx_global + 2 * nx*ny;

						//zentraler Differezenquotient in x-Richtung möglich
						wert[0] = -2.0;
						wert[1] = 1.0;
						wert[2] = 1.0;
						//modifizierter Differenzenquotient in y/z-Richtung
						wert[0] += 2.0 * 11.0 / 38.0;
						wert[3] = -28.0 / 38.0;
						wert[4] = 17.0 / 38.0;
						wert[5] = -28.0 / 38.0;
						wert[6] = 17.0 / 38.0;

						////zentraler Differezenquotient in x-Richtung möglich
						//zeile[vtx_global] = -2.0;
						//zeile[vtx_global+1] = 1.0;
						//zeile[vtx_global-1] = 1.0;
						////modifizierter Differenzenquotient in y/z-Richtung
						//zeile[vtx_global] += 2.0*11.0/38.0;
						//zeile[vtx_global+nx] = -28.0/38.0;
						//zeile[vtx_global+2*nx] = 17.0/38.0;
						//zeile[vtx_global+nx*ny] = -28.0/38.0;
						//zeile[vtx_global+2*nx*ny] = 17.0/38.0;

						//NeumannRB, Normalenvektor ist (0,1/sqrt(2),1/sqrt(2))
						//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten

						wert[0] += 2.0 * 1.0 / sqrt(2.0) * 3.0 / 2.0 * h;
						wert[3] += 1.0 / sqrt(2.0)*(-h) / 2.0;
						wert[4] += 1.0 / sqrt(2.0) * 2.0 * h;
						wert[5] += 1.0 / sqrt(2.0)*(-h) / 2.0;
						wert[6] += 1.0 / sqrt(2.0) * 2.0 * h;

						//zeile[vtx_global] += 2.0*1.0/sqrt(2.0)*3.0/2.0*h;
						//zeile[vtx_global+nx] += 1.0/sqrt(2.0)*(-h)/2.0;
						//zeile[vtx_global+2*nx] += 1.0/sqrt(2.0) *2.0*h;
						//zeile[vtx_global+nx*ny] += 1.0/sqrt(2.0)*(-h)/2.0;
						//zeile[vtx_global+2*nx*ny] += 1.0/sqrt(2.0) *2.0*h;
					}

					else if (i > nx*(ny - 1)) //hintere untere Kante ohne Eckpunkt
					{
						index[0] = vtx_global;
						index[1] = vtx_global + 1;
						index[2] = vtx_global - 1;
						index[3] = vtx_global - nx;
						index[4] = vtx_global - 2 * nx;
						index[5] = vtx_global + nx*ny;
						index[6] = vtx_global + 2 * nx*ny;

						//zentraler Differezenquotient in x-Richtung möglich
						wert[0] = -2.0;
						wert[1] = 1.0;
						wert[2] = 1.0;
						//modifizierter Differenzenquotient in y/z-Richtung
						wert[0] += 2.0 * 11.0 / 38.0;
						wert[3] = -28.0 / 38.0;
						wert[4] = 17.0 / 38.0;
						wert[5] = -28.0 / 38.0;
						wert[6] = 17.0 / 38.0;

						////zentraler Differezenquotient in x-Richtung möglich
						//zeile[vtx_global] = -2.0;
						//zeile[vtx_global+1] = 1.0;
						//zeile[vtx_global-1] = 1.0;
						//modifizierter Differenzenquotient in y/z-Richtung
						//zeile[vtx_global] += 2.0*11.0/38.0;
						//zeile[vtx_global-nx] = -28.0/38.0;
						//zeile[vtx_global-2*nx] = 17.0/38.0;
						//zeile[vtx_global+nx*ny] = -28.0/38.0;
						//zeile[vtx_global+2*nx*ny] = 17.0/38.0;

						//NeumannRB, Normalenvektor ist (0,-1/sqrt(2),1/sqrt(2))
						//RB wird auf die normale Zeile addiert, um die quadratische Struktur 
						//beizubehalten

						//wert[0] += 0.0*1.0/sqrt(2.0)*3.0/2.0*h;
						wert[3] += (-1.0) * 1.0 / sqrt(2.0)*(-h) / 2.0;
						wert[4] += (-1.0) * 1.0 / sqrt(2.0) * 2.0 * h;
						wert[5] += 1.0 / sqrt(2.0)*(-h) / 2.0;
						wert[6] += 1.0 / sqrt(2.0) * 2.0 * h;

						//zeile[vtx_global] += 0*1.0/sqrt(2.0)*3.0/2.0*h;
						//zeile[vtx_global-nx] += (-1.0)*1.0/sqrt(2.0)*(-h)/2.0;
						//zeile[vtx_global-2*nx] += (-1.0)*1.0/sqrt(2.0) *2.0*h;
						//zeile[vtx_global+nx*ny] += 1.0/sqrt(2.0)*(-h)/2.0;
						//zeile[vtx_global+2*nx*ny] += 1.0/sqrt(2.0) *2.0*h;
					}

					else if (i % nx == 1) //linke untere Kante ohne Eckpunkt
					{
						index[0] = vtx_global;
						index[1] = vtx_global + nx;
						index[2] = vtx_global - nx;
						index[3] = vtx_global + 1;
						index[4] = vtx_global + 2;
						index[5] = vtx_global + nx*ny;
						index[6] = vtx_global + 2 * nx*ny;

						//zentraler Differezenquotient in y-Richtung möglich
						wert[0] = -2.0;
						wert[1] = 1.0;
						wert[2] = 1.0;
						//modifizierter Differenzenquotient in x/z-Richtung
						wert[0] += 2.0 * 11.0 / 38.0;
						wert[3] = -28.0 / 38.0;
						wert[4] = 17.0 / 38.0;
						wert[5] = -28.0 / 38.0;
						wert[6] = 17.0 / 38.0;

						////zentraler Differezenquotient in y-Richtung möglich
						//zeile[vtx_global] = -2.0;
						//zeile[vtx_global+nx] = 1.0;
						//zeile[vtx_global-nx] = 1.0;
						////modifizierter Differenzenquotient in x/z-Richtung
						//zeile[vtx_global] += 2.0*11.0/38.0;
						//zeile[vtx_global+1] = -28.0/38.0;
						//zeile[vtx_global+2] = 17.0/38.0;
						//zeile[vtx_global+nx*ny] = -28.0/38.0;
						//zeile[vtx_global+2*nx*ny] = 17.0/38.0;

						//NeumannRB, Normalenvektor ist (1/sqrt(2),0,1/sqrt(2))
						//RB wird auf die normale Zeile addiert, um die quadratische Struktur 
						//beizubehalten

						wert[0] += 2.0 * 1.0 / sqrt(2.0) * 3.0 / 2.0 * h;
						wert[3] += 1.0 / sqrt(2.0)*(-h) / 2.0;
						wert[4] += 1.0 / sqrt(2.0) * 2.0 * h;
						wert[5] += 1.0 / sqrt(2.0)*(-h) / 2.0;
						wert[6] += 1.0 / sqrt(2.0) * 2.0 * h;

						//zeile[vtx_global] += 2.0*1.0/sqrt(2.0)*3.0/2.0*h;
						//zeile[vtx_global+1] += 1.0/sqrt(2.0)*(-h)/2.0;
						//zeile[vtx_global+2] += 1.0/sqrt(2.0) *2.0*h;
						//zeile[vtx_global+nx*ny] += 1.0/sqrt(2.0)*(-h)/2.0;
						//zeile[vtx_global+2*nx*ny] += 1.0/sqrt(2.0) *2.0*h;
					}

					else if (i % nx == 0) //rechte untere Kante ohne Eckpunkt
					{
						index[0] = vtx_global;
						index[1] = vtx_global + nx;
						index[2] = vtx_global - nx;
						index[3] = vtx_global - 1;
						index[4] = vtx_global - 2;
						index[5] = vtx_global + nx*ny;
						index[6] = vtx_global + 2 * nx*ny;

						//zentraler Differezenquotient in y-Richtung möglich
						wert[0] = -2.0;
						wert[1] = 1.0;
						wert[2] = 1.0;
						//modifizierter Differenzenquotient in x/z-Richtung
						wert[0] += 2.0 * 11.0 / 38.0;
						wert[3] = -28.0 / 38.0;
						wert[4] = 17.0 / 38.0;
						wert[5] = -28.0 / 38.0;
						wert[6] = 17.0 / 38.0;

						////zentraler Differezenquotient in y-Richtung möglich
						//zeile[vtx_global] = -2.0;
						//zeile[vtx_global+nx] = 1.0;
						//zeile[vtx_global-nx] = 1.0;
						////modifizierter Differenzenquotient in x/z-Richtung
						//zeile[vtx_global] += 2.0*11.0/38.0;
						//zeile[vtx_global-1] = -28.0/38.0;
						//zeile[vtx_global-2] = 17.0/38.0;
						//zeile[vtx_global+nx*ny] = -28.0/38.0;
						//zeile[vtx_global+2*nx*ny] = 17.0/38.0;

						//NeumannRB, Normalenvektor ist (-1/sqrt(2),0,1/sqrt(2))
						//RB wird auf die normale Zeile addiert, um die quadratische Struktur 
						//beizubehalten

						//wert[0] += 0.0*1.0/sqrt(2.0)*3.0/2.0*h;
						wert[3] += (-1.0) * 1.0 / sqrt(2.0)*(-h) / 2.0;
						wert[4] += (-1.0) * 1.0 / sqrt(2.0) * 2.0 * h;
						wert[5] += 1.0 / sqrt(2.0)*(-h) / 2.0;
						wert[6] += 1.0 / sqrt(2.0) * 2.0 * h;

						////zeile[vtx_global] += 0*1.0/sqrt(2.0)*3.0/2.0*h;
						//zeile[vtx_global-1] += (-1.0)*1.0/sqrt(2.0)*(-h)/2.0;
						//zeile[vtx_global-2] += (-1.0)*1.0/sqrt(2.0) *2.0*h;
						//zeile[vtx_global+nx*ny] += 1.0/sqrt(2.0)*(-h)/2.0;
						//zeile[vtx_global+2*nx*ny] += 1.0/sqrt(2.0) *2.0*h;
					}

					else // "innere" Punkte des Bodens
					{
						index[0] = vtx_global;
						index[1] = vtx_global + 1;
						index[2] = vtx_global - 1;
						index[3] = vtx_global + nx;
						index[4] = vtx_global - nx;
						index[5] = vtx_global + nx*ny;
						index[6] = vtx_global + 2 * nx*ny;

						////zentraler Differenzenquotient in x/y-Richtung möglich

						wert[0] = -4.0;
						wert[1] = 1.0;
						wert[2] = 1.0;
						wert[3] = 1.0;
						wert[4] = 1.0;
						//modifizierter Differenzenquotient in z-Richtung
						wert[0] += 11.0 / 38.0;
						wert[5] = -28.0 / 38.0;
						wert[6] = 17.0 / 38.0;

						////zentraler Differenzenquotient in x/y-Richtung möglich
						//zeile[vtx_global] = -4.0;
						//zeile[vtx_global+1] = 1.0;
						//zeile[vtx_global-1] = 1.0;
						//zeile[vtx_global+nx] = 1.0;
						//zeile[vtx_global-nx] = 1.0;
						//modifizierter Differenzenquotient in z-Richtung
						//zeile[vtx_global] += 11.0/38.0;
						//zeile[vtx_global+nx*ny] = -28.0/38.0;
						//zeile[vtx_global+2*nx*ny] = 17.0/38.0;

						//NeumannRB, Normalenvektor ist (0,0,1))
						//RB wird auf die normale Zeile addiert, um die quadratische Struktur 
						//beizubehalten

						wert[0] += 3.0 / 2.0 * h;
						wert[5] += (-h) / 2.0;
						wert[6] += 2.0 * h;

						//zeile[vtx_global] += 3.0/2.0*h;
						//zeile[vtx_global+nx*ny] += (-h)/2.0;
						//zeile[vtx_global+2*nx*ny] += 2.0*h;

					}
				}

				else if (i > nx*ny*(nz - 1)) //Deckel
				{
					if (i == nx*ny*(nz - 1)+1) //vorderer oberer linker Eckpunkt
					{
						index[0] = vtx_global;
						index[1] = vtx_global + 1;
						index[2] = vtx_global + 2;
						index[3] = vtx_global + nx;
						index[4] = vtx_global + 2 * nx;
						index[5] = vtx_global - nx*ny;
						index[6] = vtx_global - 2 * nx*ny;

						//zentraler Differenzenquotient in keine Richtung möglich
						wert[0] = 3.0 * 11.0 / 38.0;
						wert[1] = -28.0 / 38.0;
						wert[2] = 17.0 / 38.0;
						wert[3] = -28.0 / 38.0;
						wert[4] = 17.0 / 38.0;
						wert[5] = -28.0 / 38.0;
						wert[6] = 17.0 / 38.0;

						////zentraler Differenzenquotient in keine Richtung möglich
						//zeile[vtx_global] = 3.0*11.0/38.0;
						//zeile[vtx_global+1] = -28.0/38.0;
						//zeile[vtx_global+2] = 17.0/38.0;
						//zeile[vtx_global+nx] = -28.0/38.0;
						//zeile[vtx_global+2*nx] = 17.0/38.0;
						//zeile[vtx_global-nx*ny] = -28.0/38.0;
						//zeile[vtx_global-2*nx*ny] = 17.0/38.0;

						//NeumannRB, Normalenvektor ist (1/sqrt(3),1/sqrt(3),-1/sqrt(3))
						//RB wird auf die normale Zeile addiert, um die quadratische Struktur 
						//beizubehalten

						wert[0] += 1.0 * 1.0 / sqrt(3.0) * 3.0 / 2.0 * h;
						wert[1] += 1.0 / sqrt(3.0)*(-h) / 2.0;
						wert[2] += 1.0 / sqrt(3.0) * 2.0 * h;
						wert[3] += 1.0 / sqrt(3.0)*(-h) / 2.0;
						wert[4] += 1.0 / sqrt(3.0) * 2.0 * h;
						wert[5] += (-1.0) * 1.0 / sqrt(3.0)*(-h) / 2.0;
						wert[6] += (-1.0) * 1.0 / sqrt(3.0) * 2.0 * h;

						//zeile[vtx_global] += 1.0*1.0/sqrt(3.0)*3.0/2.0*h;
						//zeile[vtx_global+1] += 1.0/sqrt(3.0)*(-h)/2.0;
						//zeile[vtx_global+2] += 1.0/sqrt(3.0) *2.0*h;
						//zeile[vtx_global+nx] += 1.0/sqrt(3.0)*(-h)/2.0;
						//zeile[vtx_global+2*nx] += 1.0/sqrt(3.0) *2.0*h;
						//zeile[vtx_global-nx*ny] += (-1.0)*1.0/sqrt(3.0)*(-h)/2.0;
						//zeile[vtx_global-2*nx*ny] += (-1.0)*1.0/sqrt(3.0) *2.0*h;
					}

					else if (i == nx*ny*(nz - 1)+nx) //vorderer oberer rechter Eckpunkt
					{
						index[0] = vtx_global;
						index[1] = vtx_global - 1;
						index[2] = vtx_global - 2;
						index[3] = vtx_global + nx;
						index[4] = vtx_global + 2 * nx;
						index[5] = vtx_global - nx*ny;
						index[6] = vtx_global - 2 * nx*ny;

						//zentraler Differenzenquotient in keine Richtung möglich
						wert[0] = 3.0 * 11.0 / 38.0;
						wert[1] = -28.0 / 38.0;
						wert[2] = 17.0 / 38.0;
						wert[3] = -28.0 / 38.0;
						wert[4] = 17.0 / 38.0;
						wert[5] = -28.0 / 38.0;
						wert[6] = 17.0 / 38.0;

						////zentraler Differenzenquotient in keine Richtung möglich
						//zeile[vtx_global] = 3.0*11.0/38.0;
						//zeile[vtx_global-1] = -28.0/38.0;
						//zeile[vtx_global-2] = 17.0/38.0;
						//zeile[vtx_global+nx] = -28.0/38.0;
						//zeile[vtx_global+2*nx] = 17.0/38.0;
						//zeile[vtx_global-nx*ny] = -28.0/38.0;
						//zeile[vtx_global-2*nx*ny] = 17.0/38.0;

						//NeumannRB, Normalenvektor ist (-1/sqrt(3),1/sqrt(3),-1/sqrt(3))
						//RB wird auf die normale Zeile addiert, um die quadratische Struktur 
						//beizubehalten

						wert[0] += (-1.0) * 1.0 / sqrt(3.0) * 3.0 / 2.0 * h;
						wert[1] += (-1.0) * 1.0 / sqrt(3.0)*(-h) / 2.0;
						wert[2] += (-1.0) * 1.0 / sqrt(3.0) * 2.0 * h;
						wert[3] += 1.0 / sqrt(3.0)*(-h) / 2.0;
						wert[4] += 1.0 / sqrt(3.0) * 2.0 * h;
						wert[5] += (-1.0) * 1.0 / sqrt(3.0)*(-h) / 2.0;
						wert[6] += (-1.0) * 1.0 / sqrt(3.0) * 2.0 * h;

						//zeile[vtx_global] += (-1.0)*1.0/sqrt(3.0)*3.0/2.0*h;
						//zeile[vtx_global-1] += (-1.0)*1.0/sqrt(3.0)*(-h)/2.0;
						//zeile[vtx_global-2] += (-1.0)*1.0/sqrt(3.0) *2.0*h;
						//zeile[vtx_global+nx] += 1.0/sqrt(3.0)*(-h)/2.0;
						//zeile[vtx_global+2*nx] += 1.0/sqrt(3.0) *2.0*h;
						//zeile[vtx_global-nx*ny] += (-1.0)*1.0/sqrt(3.0)*(-h)/2.0;
						//zeile[vtx_global-2*nx*ny] += (-1.0)*1.0/sqrt(3.0) *2.0*h;
					}

					else if (i == nx*ny*nz-nx+1) //hinterer oberer linker Eckpunkt
					{
						index[0] = vtx_global;
						index[1] = vtx_global + 1;
						index[2] = vtx_global + 2;
						index[3] = vtx_global - nx;
						index[4] = vtx_global - 2 * nx;
						index[5] = vtx_global - nx*ny;
						index[6] = vtx_global - 2 * nx*ny;

						//zentraler Differenzenquotient in keine Richtung möglich
						wert[0] = 3.0 * 11.0 / 38.0;
						wert[1] = -28.0 / 38.0;
						wert[2] = 17.0 / 38.0;
						wert[3] = -28.0 / 38.0;
						wert[4] = 17.0 / 38.0;
						wert[5] = -28.0 / 38.0;
						wert[6] = 17.0 / 38.0;

						////zentraler Differenzenquotient in keine Richtung möglich
						//zeile[vtx_global] = 3.0*11.0/38.0;
						//zeile[vtx_global+1] = -28.0/38.0;
						//zeile[vtx_global+2] = 17.0/38.0;
						//zeile[vtx_global-nx] = -28.0/38.0;
						//zeile[vtx_global-2*nx] = 17.0/38.0;
						//zeile[vtx_global-nx*ny] = -28.0/38.0;
						//zeile[vtx_global-2*nx*ny] = 17.0/38.0;

						//NeumannRB, Normalenvektor ist (1/sqrt(3),-1/sqrt(3),-1/sqrt(3))
						//RB wird auf die normale Zeile addiert, um die quadratische Struktur 
						//beizubehalten

						wert[0] += (-1.0) * 1.0 / sqrt(3.0) * 3.0 / 2.0 * h;
						wert[1] += 1.0 / sqrt(3.0)*(-h) / 2.0;
						wert[2] += 1.0 / sqrt(3.0) * 2.0 * h;
						wert[3] += (-1.0) * 1.0 / sqrt(3.0)*(-h) / 2.0;
						wert[4] += (-1.0) * 1.0 / sqrt(3.0) * 2.0 * h;
						wert[5] += (-1.0) * 1.0 / sqrt(3.0)*(-h) / 2.0;
						wert[6] += (-1.0) * 1.0 / sqrt(3.0) * 2.0 * h;

						//zeile[vtx_global] += (-1.0)*1.0/sqrt(3.0)*3.0/2.0*h;
						//zeile[vtx_global+1] += 1.0/sqrt(3.0)*(-h)/2.0;
						//zeile[vtx_global+2] += 1.0/sqrt(3.0) *2.0*h;
						//zeile[vtx_global-nx] += (-1.0)*1.0/sqrt(3.0)*(-h)/2.0;
						//zeile[vtx_global-2*nx] += (-1.0)*1.0/sqrt(3.0) *2.0*h;
						//zeile[vtx_global-nx*ny] += (-1.0)*1.0/sqrt(3.0)*(-h)/2.0;
						//zeile[vtx_global-2*nx*ny] += (-1.0)*1.0/sqrt(3.0) *2.0*h;
					}

					else if (i == nx*ny*nz) //hinterer  oberer rechter Eckpunkt
					{
						index[0] = vtx_global;
						index[1] = vtx_global - 1;
						index[2] = vtx_global - 2;
						index[3] = vtx_global - nx;
						index[4] = vtx_global - 2 * nx;
						index[5] = vtx_global - nx*ny;
						index[6] = vtx_global - 2 * nx*ny;

						//zentraler Differenzenquotient in keine Richtung möglich
						wert[0] = 3.0 * 11.0 / 38.0;
						wert[1] = -28.0 / 38.0;
						wert[2] = 17.0 / 38.0;
						wert[3] = -28.0 / 38.0;
						wert[4] = 17.0 / 38.0;
						wert[5] = -28.0 / 38.0;
						wert[6] = 17.0 / 38.0;

						////zentraler Differenzenquotient in keine Richtung möglich
						//zeile[vtx_global] = 3.0*11.0/38.0;
						//zeile[vtx_global-1] = -28.0/38.0;
						//zeile[vtx_global-2] = 17.0/38.0;
						//zeile[vtx_global-nx] = -28.0/38.0;
						//zeile[vtx_global-2*nx] = 17.0/38.0;
						//zeile[vtx_global-nx*ny] = -28.0/38.0;
						//zeile[vtx_global-2*nx*ny] = 17.0/38.0;

						//NeumannRB, Normalenvektor ist (-1/sqrt(3),-1/sqrt(3),-1/sqrt(3))
						//RB wird auf die normale Zeile addiert, um die quadratische Struktur 
						//beizubehalten

						wert[0] += (-3.0) * 1.0 / sqrt(3.0) * 3.0 / 2.0 * h;
						wert[1] += (-1.0) * 1.0 / sqrt(3.0)*(-h) / 2.0;
						wert[2] += (-1.0) * 1.0 / sqrt(3.0) * 2.0 * h;
						wert[3] += (-1.0) * 1.0 / sqrt(3.0)*(-h) / 2.0;
						wert[4] += (-1.0) * 1.0 / sqrt(3.0) * 2.0 * h;
						wert[5] += (-1.0) * 1.0 / sqrt(3.0)*(-h) / 2.0;
						wert[6] += (-1.0) * 1.0 / sqrt(3.0) * 2.0 * h;

						//zeile[vtx_global] += (-3.0)*1.0/sqrt(3.0)*3.0/2.0*h;
						//zeile[vtx_global-1] += (-1.0)*1.0/sqrt(3.0)*(-h)/2.0;
						//zeile[vtx_global-2] += (-1.0)*1.0/sqrt(3.0) *2.0*h;
						//zeile[vtx_global-nx] += (-1.0)*1.0/sqrt(3.0)*(-h)/2.0;
						//zeile[vtx_global-2*nx] += (-1.0)*1.0/sqrt(3.0) *2.0*h;
						//zeile[vtx_global-nx*ny] += (-1.0)*1.0/sqrt(3.0)*(-h)/2.0;
						//zeile[vtx_global-2*nx*ny] += (-1.0)*1.0/sqrt(3.0) *2.0*h;
					}

					else if (i  < nx*ny*(nz - 1)+nx) //vordere obere Kante ohne Eckpunkt
					{
						index[0] = vtx_global;
						index[1] = vtx_global + 1;
						index[2] = vtx_global - 1;
						index[3] = vtx_global + nx;
						index[4] = vtx_global + 2 * nx;
						index[5] = vtx_global - nx*ny;
						index[6] = vtx_global - 2 * nx*ny;

						//zentraler Differezenquotient in x-Richtung möglich
						wert[0] = -2.0;
						wert[1] = 1.0;
						wert[2] = 1.0;
						//modifizierter Differenzenquotient in y/z-Richtung
						wert[0] += 2.0 * 11.0 / 38.0;
						wert[3] = -28.0 / 38.0;
						wert[4] = 17.0 / 38.0;
						wert[5] = -28.0 / 38.0;
						wert[6] = 17.0 / 38.0;

						////zentraler Differezenquotient in x-Richtung möglich
						//zeile[vtx_global] = -2.0;
						//zeile[vtx_global+1] = 1.0;
						//zeile[vtx_global-1] = 1.0;
						////modifizierter Differenzenquotient in y/z-Richtung
						//zeile[vtx_global] += 2.0*11.0/38.0;
						//zeile[vtx_global+nx] = -28.0/38.0;
						//zeile[vtx_global+2*nx] = 17.0/38.0;
						//zeile[vtx_global-nx*ny] = -28.0/38.0;
						//zeile[vtx_global-2*nx*ny] = 17.0/38.0;

						//NeumannRB, Normalenvektor ist (0,1/sqrt(2),-1/sqrt(2))
						//RB wird auf die normale Zeile addiert, um die quadratische Struktur 
						//beizubehalten

						//wert[0] += 0.0*1.0/sqrt(2.0)*3.0/2.0*h;
						wert[3] += 1.0 / sqrt(2.0)*(-h) / 2.0;
						wert[4] += 1.0 / sqrt(2.0) * 2.0 * h;
						wert[5] += (-1.0) * 1.0 / sqrt(2.0)*(-h) / 2.0;
						wert[6] += (-1.0) * 1.0 / sqrt(2.0) * 2.0 * h;

						////zeile[vtx_global] += 0.0*1.0/sqrt(2.0)*3.0/2.0*h;
						//zeile[vtx_global+nx] += 1.0/sqrt(2.0)*(-h)/2.0;
						//zeile[vtx_global+2*nx] += 1.0/sqrt(2.0) *2.0*h;
						//zeile[vtx_global-nx*ny] += (-1.0)*1.0/sqrt(2.0)*(-h)/2.0;
						//zeile[vtx_global-2*nx*ny] += (-1.0)*1.0/sqrt(2.0) *2.0*h;
					}

					else if (i > nx*ny*(nz - 1)+nx*(ny-1)) //hintere obere Kante ohne Eckpunkt
					{
						index[0] = vtx_global;
						index[1] = vtx_global + 1;
						index[2] = vtx_global - 1;
						index[3] = vtx_global - nx;
						index[4] = vtx_global - 2 * nx;
						index[5] = vtx_global - nx*ny;
						index[6] = vtx_global - 2 * nx*ny;

						//zentraler Differezenquotient in x-Richtung möglich
						wert[0] = -2.0;
						wert[1] = 1.0;
						wert[2] = 1.0;
						//modifizierter Differenzenquotient in y/z-Richtung
						wert[0] += 2.0 * 11.0 / 38.0;
						wert[3] = -28.0 / 38.0;
						wert[4] = 17.0 / 38.0;
						wert[5] = -28.0 / 38.0;
						wert[6] = 17.0 / 38.0;

						////zentraler Differezenquotient in x-Richtung möglich
						//zeile[vtx_global] = -2.0;
						//zeile[vtx_global+1] = 1.0;
						//zeile[vtx_global-1] = 1.0;
						////modifizierter Differenzenquotient in y/z-Richtung
						//zeile[vtx_global] += 2.0*11.0/38.0;
						//zeile[vtx_global-nx] = -28.0/38.0;
						//zeile[vtx_global-2*nx] = 17.0/38.0;
						//zeile[vtx_global-nx*ny] = -28.0/38.0;
						//zeile[vtx_global-2*nx*ny] = 17.0/38.0;

						//NeumannRB, Normalenvektor ist (0,-1/sqrt(2),-1/sqrt(2))
						//RB wird auf die normale Zeile addiert, um die quadratische Struktur 
						//beizubehalten

						wert[0] += (-2.0) * 1.0 / sqrt(2.0) * 3.0 / 2.0 * h;
						wert[3] += (-1.0) * 1.0 / sqrt(2.0)*(-h) / 2.0;
						wert[4] += (-1.0) * 1.0 / sqrt(2.0) * 2.0 * h;
						wert[5] += (-1.0) * 1.0 / sqrt(2.0)*(-h) / 2.0;
						wert[6] += (-1.0) * 1.0 / sqrt(2.0) * 2.0 * h;

						//zeile[vtx_global] += (-2.0)*1.0/sqrt(2.0)*3.0/2.0*h;
						//zeile[vtx_global-nx] += (-1.0)*1.0/sqrt(2.0)*(-h)/2.0;
						//zeile[vtx_global-2*nx] += (-1.0)*1.0/sqrt(2.0) *2.0*h;
						//zeile[vtx_global-nx*ny] += (-1.0)*1.0/sqrt(2.0)*(-h)/2.0;
						//zeile[vtx_global-2*nx*ny] += (-1.0)*1.0/sqrt(2.0) *2.0*h;
					}

					else if (i % nx== 1) //linke obere Kante ohne Eckpunkt
					{
						index[0] = vtx_global;
						index[1] = vtx_global + nx;
						index[2] = vtx_global - nx;
						index[3] = vtx_global + 1;
						index[4] = vtx_global + 2;
						index[5] = vtx_global - nx*ny;
						index[6] = vtx_global - 2 * nx*ny;

						//zentraler Differezenquotient in y-Richtung möglich
						wert[0] = -2.0;
						wert[1] = 1.0;
						wert[2] = 1.0;
						//modifizierter Differenzenquotient in x/z-Richtung
						wert[0] += 2.0 * 11.0 / 38.0;
						wert[3] = -28.0 / 38.0;
						wert[4] = 17.0 / 38.0;
						wert[5] = -28.0 / 38.0;
						wert[6] = 17.0 / 38.0;

						////zentraler Differezenquotient in y-Richtung möglich
						//zeile[vtx_global] = -2.0;
						//zeile[vtx_global+nx] = 1.0;
						//zeile[vtx_global-nx] = 1.0;
						////modifizierter Differenzenquotient in x/z-Richtung
						//zeile[vtx_global] += 2.0*11.0/38.0;
						//zeile[vtx_global+1] = -28.0/38.0;
						//zeile[vtx_global+2] = 17.0/38.0;
						//zeile[vtx_global-nx*ny] = -28.0/38.0;
						//zeile[vtx_global-2*nx*ny] = 17.0/38.0;

						//NeumannRB, Normalenvektor ist (1/sqrt(2),0,-1/sqrt(2))
						//RB wird auf die normale Zeile addiert, um die quadratische Struktur 
						//beizubehalten

						//wert[0] += 0.0*1.0/sqrt(2.0)*3.0/2.0*h;
						wert[3] += 1.0 / sqrt(2.0)*(-h) / 2.0;
						wert[4] += 1.0 / sqrt(2.0) * 2.0 * h;
						wert[5] += (-1.0) * 1.0 / sqrt(2.0)*(-h) / 2.0;
						wert[6] += (-1.0) * 1.0 / sqrt(2.0) * 2.0 * h;

						////zeile[vtx_global] += 0*1.0/sqrt(2.0)*3.0/2.0*h;
						//zeile[vtx_global+1] += 1.0/sqrt(2.0)*(-h)/2.0;
						//zeile[vtx_global+2] += 1.0/sqrt(2.0) *2.0*h;
						//zeile[vtx_global-nx*ny] += (-1.0)*1.0/sqrt(2.0)*(-h)/2.0;
						//zeile[vtx_global-2*nx*ny] += (-1.0)*1.0/sqrt(2.0) *2.0*h;
					}

					else if (i % nx == 0) //rechte obere Kante ohne Eckpunkt
					{
						index[0] = vtx_global;
						index[1] = vtx_global + nx;
						index[2] = vtx_global - nx;
						index[3] = vtx_global - 1;
						index[4] = vtx_global - 2;
						index[5] = vtx_global - nx*ny;
						index[6] = vtx_global - 2 * nx*ny;

						//zentraler Differezenquotient in y-Richtung möglich
						wert[0] = -2.0;
						wert[1] = 1.0;
						wert[2] = 1.0;
						//modifizierter Differenzenquotient in x/z-Richtung
						wert[0] += 2.0 * 11.0 / 38.0;
						wert[3] = -28.0 / 38.0;
						wert[4] = 17.0 / 38.0;
						wert[5] = -28.0 / 38.0;
						wert[6] = 17.0 / 38.0;

						////zentraler Differezenquotient in y-Richtung möglich
						//zeile[vtx_global] = -2.0;
						//zeile[vtx_global+nx] = 1.0;
						//zeile[vtx_global-nx] = 1.0;
						////modifizierter Differenzenquotient in x/z-Richtung
						//zeile[vtx_global] += 2.0*11.0/38.0;
						//zeile[vtx_global-1] = -28.0/38.0;
						//zeile[vtx_global-2] = 17.0/38.0;
						//zeile[vtx_global-nx*ny] = -28.0/38.0;
						//zeile[vtx_global-2*nx*ny] = 17.0/38.0;

						//NeumannRB, Normalenvektor ist (-1/sqrt(2),0,-1/sqrt(2))
						//RB wird auf die normale Zeile addiert, um die quadratische Struktur 
						//beizubehalten

						wert[0] += (-2.0) * 1.0 / sqrt(2.0) * 3.0 / 2.0 * h;
						wert[3] += (-1.0) * 1.0 / sqrt(2.0)*(-h) / 2.0;
						wert[4] += (-1.0) * 1.0 / sqrt(2.0) * 2.0 * h;
						wert[5] += (-1.0) * 1.0 / sqrt(2.0)*(-h) / 2.0;
						wert[6] += (-1.0) * 1.0 / sqrt(2.0) * 2.0 * h;

						//zeile[vtx_global] += (-2.0)*1.0/sqrt(2.0)*3.0/2.0*h;
						//zeile[vtx_global-1] += (-1.0)*1.0/sqrt(2.0)*(-h)/2.0;
						//zeile[vtx_global-2] += (-1.0)*1.0/sqrt(2.0) *2.0*h;
						//zeile[vtx_global-nx*ny] += (-1.0)*1.0/sqrt(2.0)*(-h)/2.0;
						//zeile[vtx_global-2*nx*ny] += (-1.0)*1.0/sqrt(2.0) *2.0*h;
					}

					else // "innere" Punkte des Deckels
					{
						index[0] = vtx_global;
						index[1] = vtx_global + 1;
						index[2] = vtx_global - 1;
						index[3] = vtx_global + nx;
						index[4] = vtx_global - nx;
						index[5] = vtx_global - nx*ny;
						index[6] = vtx_global - 2 * nx*ny;

						//zentraler Differenzenquotient in x/y-Richtung möglich
						wert[0] = -4.0;
						wert[1] = 1.0;
						wert[2] = 1.0;
						wert[3] = 1.0;
						wert[4] = 1.0;
						//modifizierter Differenzenquotient in z-Richtung
						wert[0] += 11.0 / 38.0;
						wert[5] = -28.0 / 38.0;
						wert[6] = 17.0 / 38.0;

						////zentraler Differenzenquotient in x/y-Richtung möglich
						//zeile[vtx_global] = -4.0;
						//zeile[vtx_global+1] = 1.0;
						//zeile[vtx_global-1] = 1.0;
						//zeile[vtx_global+nx] = 1.0;
						//zeile[vtx_global-nx] = 1.0;
						////modifizierter Differenzenquotient in z-Richtung
						//zeile[vtx_global] += 11.0/38.0;
						//zeile[vtx_global-nx*ny] = -28.0/38.0;
						//zeile[vtx_global-2*nx*ny] = 17.0/38.0;

						//NeumannRB, Normalenvektor ist (0,0,-1))
						//RB wird auf die normale Zeile addiert, um die quadratische Struktur 
						//beizubehalten

						wert[0] += (-1.0) * 3.0 / 2.0 * h;
						wert[5] += (-1.0)*(-h) / 2.0;
						wert[6] += (-1.0) * 2.0 * h;

						//zeile[vtx_global] += (-1.0)*3.0/2.0*h;
						//zeile[vtx_global-nx*ny] += (-1.0)*(-h)/2.0;
						//zeile[vtx_global-2*nx*ny] += (-1.0)*2.0*h;

					}
				}

				else if (i % (nx*ny) <= nx) //vordere Seite, aber nicht Boden oder Deckel
				{
					if (i % nx == 1) //linke Kante
					{
						index[0] = vtx_global;
						index[1] = vtx_global + nx*ny;
						index[2] = vtx_global - nx*ny;
						index[3] = vtx_global + 1;
						index[4] = vtx_global + 2;
						index[5] = vtx_global + nx;
						index[6] = vtx_global + 2 * nx;

						//zentraler Differenzenquotient nur in z-Richtung möglich
						wert[0] = -2.0;
						wert[1] = 1.0;
						wert[2] = 1.0;
						//modifizierter Differenzenquotient in x/y-Richtung
						wert[0] += 2.0 * 11.0 / 38.0;
						wert[3] = -28.0 / 38.0;
						wert[4] = 17.0 / 38.0;
						wert[5] = -28.0 / 38.0;
						wert[6] = 17.0 / 38.0;

						////zentraler Differenzenquotient nur in z-Richtung möglich
						//zeile[vtx_global] = -2.0;
						//zeile[vtx_global+nx*ny] = 1.0;
						//zeile[vtx_global-nx*ny] = 1.0;
						////modifizierter Differenzenquotient in x/y-Richtung
						//zeile[vtx_global] += 2.0*11.0/38.0;
						//zeile[vtx_global+1] = -28.0/38.0;
						//zeile[vtx_global+2] = 17.0/38.0;
						//zeile[vtx_global+nx] = -28.0/38.0;
						//zeile[vtx_global+2*nx] = 17.0/38.0;

						//NeumannRB, Normalenvektor ist (1/sqrt(2),1/sqrt(2),0)
						//RB wird auf die normale Zeile addiert, um die quadratische Struktur 
						//beizubehalten

						wert[0] += 2.0 * 1.0 / sqrt(2.0) * 3.0 / 2.0 * h;
						wert[3] += 1.0 / sqrt(2.0)*(-h) / 2.0;
						wert[4] += 1.0 / sqrt(2.0) * 2.0 * h;
						wert[5] += 1.0 / sqrt(2.0)*(-h) / 2.0;
						wert[6] += 1.0 / sqrt(2.0) * 2.0 * h;

						//zeile[vtx_global] += 2.0*1.0/sqrt(2.0)*3.0/2.0*h;
						//zeile[vtx_global+1] += 1.0/sqrt(2.0)*(-h)/2.0;
						//zeile[vtx_global+2] += 1.0/sqrt(2.0) *2.0*h;
						//zeile[vtx_global+nx] += 1.0/sqrt(2.0)*(-h)/2.0;
						//zeile[vtx_global+2*nx] += 1.0/sqrt(2.0) *2.0*h;

					}

					if (i % nx == 0) //rechte Kante
					{
						index[0] = vtx_global;
						index[1] = vtx_global + nx*ny;
						index[2] = vtx_global - nx*ny;
						index[3] = vtx_global - 1;
						index[4] = vtx_global - 2;
						index[5] = vtx_global + nx;
						index[6] = vtx_global + 2 * nx;

						//zentraler Differenzenquotient nur in z-Richtung möglich
						wert[0] = -2.0;
						wert[1] = 1.0;
						wert[2] = 1.0;
						//modifizierter Differenzenquotient in x/y-Richtung
						wert[0] += 2.0 * 11.0 / 38.0;
						wert[3] = -28.0 / 38.0;
						wert[4] = 17.0 / 38.0;
						wert[5] = -28.0 / 38.0;
						wert[6] = 17.0 / 38.0;

						////zentraler Differenzenquotient nur in z-Richtung möglich
						//zeile[vtx_global] = -2.0;
						//zeile[vtx_global+nx*ny] = 1.0;
						//zeile[vtx_global-nx*ny] = 1.0;
						////modifizierter Differenzenquotient in x/y-Richtung
						//zeile[vtx_global] += 2.0*11.0/38.0;
						//zeile[vtx_global-1] = -28.0/38.0;
						//zeile[vtx_global-2] = 17.0/38.0;
						//zeile[vtx_global+nx] = -28.0/38.0;
						//zeile[vtx_global+2*nx] = 17.0/38.0;

						//NeumannRB, Normalenvektor ist (-1/sqrt(2),1/sqrt(2),0)
						//RB wird auf die normale Zeile addiert, um die quadratische Struktur 
						//beizubehalten

						//wert[0] += 0*1.0/sqrt(2.0)*3.0/2.0*h;
						wert[3] += (-1.0) * 1.0 / sqrt(2.0)*(-h) / 2.0;
						wert[4] += (-1.0) * 1.0 / sqrt(2.0) * 2.0 * h;
						wert[5] += 1.0 / sqrt(2.0)*(-h) / 2.0;
						wert[6] += 1.0 / sqrt(2.0) * 2.0 * h;

						////zeile[vtx_global] += 0*1.0/sqrt(2.0)*3.0/2.0*h;
						//zeile[vtx_global-1] += (-1.0)*1.0/sqrt(2.0)*(-h)/2.0;
						//zeile[vtx_global-2] += (-1.0)*1.0/sqrt(2.0) *2.0*h;
						//zeile[vtx_global+nx] += 1.0/sqrt(2.0)*(-h)/2.0;
						//zeile[vtx_global+2*nx] += 1.0/sqrt(2.0) *2.0*h;
					}

					else //vordere "innere" Seite
					{
						index[0] = vtx_global;
						index[1] = vtx_global + nx*ny;
						index[2] = vtx_global - nx*ny;
						index[3] = vtx_global + 1;
						index[4] = vtx_global - 1;
						index[5] = vtx_global + nx;
						index[6] = vtx_global + 2 * nx;

						//zentraler Differenzenquotient in x/z-Richtung möglich
						wert[0] = -4.0;
						wert[1] = 1.0;
						wert[2] = 1.0;
						wert[3] = 1.0;
						wert[4] = 1.0;
						//modifizierter Differenzenquotient in y-Richtung
						wert[0] += 11.0 / 38.0;
						wert[5] = -28.0 / 38.0;
						wert[6] = 17.0 / 38.0;

						////zentraler Differenzenquotient in x/z-Richtung möglich
						//zeile[vtx_global] = -4.0;
						//zeile[vtx_global+nx*ny] = 1.0;
						//zeile[vtx_global-nx*ny] = 1.0;
						//zeile[vtx_global+1]=1.0;
						//zeile[vtx_global-1]=1.0;
						////modifizierter Differenzenquotient in y-Richtung
						//zeile[vtx_global] += 11.0/38.0;
						//zeile[vtx_global+nx] = -28.0/38.0;
						//zeile[vtx_global+2*nx] = 17.0/38.0;

						//NeumannRB, Normalenvektor ist (0,1,0)
						//RB wird auf die normale Zeile addiert, um die quadratische Struktur 
						//beizubehalten

						wert[0] += 3.0 / 2.0 * h;
						wert[5] += (-h) / 2.0;
						wert[6] += 2.0 * h;

						//zeile[vtx_global] += 3.0/2.0*h;
						//zeile[vtx_global+nx] += (-h)/2.0;
						//zeile[vtx_global+2*nx] += 2.0*h;
					}
				}

				else if (i % (nx*ny) > nx*(ny - 1)) //hintere Seite, aber nicht Boden oder Deckel
				{
					if (i % nx == 1) //linke Kante
					{
						index[0] = vtx_global;
						index[1] = vtx_global + nx*ny;
						index[2] = vtx_global - nx*ny;
						index[3] = vtx_global + 1;
						index[4] = vtx_global + 2;
						index[5] = vtx_global - nx;
						index[6] = vtx_global - 2 * nx;

						//zentraler Differenzenquotient nur in z-Richtung möglich
						wert[0] = -2.0;
						wert[1] = 1.0;
						wert[2] = 1.0;
						//modifizierter Differenzenquotient in x/y-Richtung
						wert[0] += 2.0 * 11.0 / 38.0;
						wert[3] = -28.0 / 38.0;
						wert[4] = 17.0 / 38.0;
						wert[5] = -28.0 / 38.0;
						wert[6] = 17.0 / 38.0;

						////zentraler Differenzenquotient nur in z-Richtung möglich
						//zeile[vtx_global] = -2.0;
						//zeile[vtx_global+nx*ny] = 1.0;
						//zeile[vtx_global-nx*ny] = 1.0;
						////modifizierter Differenzenquotient in x/y-Richtung
						//zeile[vtx_global] += 2.0*11.0/38.0;
						//zeile[vtx_global+1] = -28.0/38.0;
						//zeile[vtx_global+2] = 17.0/38.0;
						//zeile[vtx_global-nx] = -28.0/38.0;
						//zeile[vtx_global-2*nx] = 17.0/38.0;

						//NeumannRB, Normalenvektor ist (1/sqrt(2),-1/sqrt(2),0)
						//RB wird auf die normale Zeile addiert, um die quadratische Struktur 
						//beizubehalten

						//wert[0] += 0.0*1.0/sqrt(2.0)*3.0/2.0*h;
						wert[3] += 1.0 / sqrt(2.0)*(-h) / 2.0;
						wert[4] += 1.0 / sqrt(2.0) * 2.0 * h;
						wert[5] += (-1.0) * 1.0 / sqrt(2.0)*(-h) / 2.0;
						wert[6] += (-1.0) * 1.0 / sqrt(2.0) * 2.0 * h;

						////zeile[vtx_global] += 0.0*1.0/sqrt(2.0)*3.0/2.0*h;
						//zeile[vtx_global+1] += 1.0/sqrt(2.0)*(-h)/2.0;
						//zeile[vtx_global+2] += 1.0/sqrt(2.0) *2.0*h;
						//zeile[vtx_global-nx] += (-1.0)*1.0/sqrt(2.0)*(-h)/2.0;
						//zeile[vtx_global-2*nx] += (-1.0)*1.0/sqrt(2.0) *2.0*h;

					}

					if (i % nx == 0) //rechte Kante
					{
						index[0] = vtx_global;
						index[1] = vtx_global + nx*ny;
						index[2] = vtx_global - nx*ny;
						index[3] = vtx_global - 1;
						index[4] = vtx_global - 2;
						index[5] = vtx_global - nx;
						index[6] = vtx_global - 2 * nx;

						//zentraler Differenzenquotient nur in z-Richtung möglich
						wert[0] = -2.0;
						wert[1] = 1.0;
						wert[2] = 1.0;
						//modifizierter Differenzenquotient in x/y-Richtung
						wert[0] += 2.0 * 11.0 / 38.0;
						wert[3] = -28.0 / 38.0;
						wert[4] = 17.0 / 38.0;
						wert[5] = -28.0 / 38.0;
						wert[6] = 17.0 / 38.0;

						////zentraler Differenzenquotient nur in z-Richtung möglich
						//zeile[vtx_global] = -2.0;
						//zeile[vtx_global+nx*ny] = 1.0;
						//zeile[vtx_global-nx*ny] = 1.0;
						////modifizierter Differenzenquotient in x/y-Richtung
						//zeile[vtx_global] += 2.0*11.0/38.0;
						//zeile[vtx_global-1] = -28.0/38.0;
						//zeile[vtx_global-2] = 17.0/38.0;
						//zeile[vtx_global-nx] = -28.0/38.0;
						//zeile[vtx_global-2*nx] = 17.0/38.0;

						//NeumannRB, Normalenvektor ist (-1/sqrt(2),-1/sqrt(2),0)
						//RB wird auf die normale Zeile addiert, um die quadratische Struktur 
						//beizubehalten
						wert[0] += (-2.0) * 1.0 / sqrt(2.0) * 3.0 / 2.0 * h;
						wert[3] += (-1.0) * 1.0 / sqrt(2.0)*(-h) / 2.0;
						wert[4] += (-1.0) * 1.0 / sqrt(2.0) * 2.0 * h;
						wert[5] += (-1.0) * 1.0 / sqrt(2.0)*(-h) / 2.0;
						wert[6] += (-1.0) * 1.0 / sqrt(2.0) * 2.0 * h;

						//zeile[vtx_global] += (-2.0)*1.0/sqrt(2.0)*3.0/2.0*h;
						//zeile[vtx_global-1] += (-1.0)*1.0/sqrt(2.0)*(-h)/2.0;
						//zeile[vtx_global-2] += (-1.0)*1.0/sqrt(2.0) *2.0*h;
						//zeile[vtx_global-nx] += (-1.0)*1.0/sqrt(2.0)*(-h)/2.0;
						//zeile[vtx_global-2*nx] += (-1.0)*1.0/sqrt(2.0) *2.0*h;
					}

					else // hintere "innere" Seite
					{
						index[0] = vtx_global;
						index[1] = vtx_global + nx*ny;
						index[2] = vtx_global - nx*ny;
						index[3] = vtx_global + 1;
						index[4] = vtx_global - 1;
						index[5] = vtx_global - nx;
						index[6] = vtx_global - 2 * nx;

						//zentraler Differenzenquotient in x/z-Richtung möglich
						wert[0] = -4.0;
						wert[1] = 1.0;
						wert[2] = 1.0;
						wert[3] = 1.0;
						wert[4] = 1.0;
						//modifizierter Differenzenquotient in y-Richtung
						wert[0] += 11.0 / 38.0;
						wert[5] = -28.0 / 38.0;
						wert[6] = 17.0 / 38.0;

						////zentraler Differenzenquotient in x/z-Richtung möglich
						//zeile[vtx_global] = -4.0;
						//zeile[vtx_global+nx*ny] = 1.0;
						//zeile[vtx_global-nx*ny] = 1.0;
						//zeile[vtx_global+1] =1.0;
						//zeile[vtx_global-1] =1.0;
						////modifizierter Differenzenquotient in y-Richtung
						//zeile[vtx_global] += 11.0/38.0;
						//zeile[vtx_global-nx] = -28.0/38.0;
						//zeile[vtx_global-2*nx] = 17.0/38.0;

						//NeumannRB, Normalenvektor ist (0,-1,0)
						//RB wird auf die normale Zeile addiert, um die quadratische Struktur 
						//beizubehalten
						wert[0] += (-1.0) * 3.0 / 2.0 * h;
						wert[5] += (-1.0)*(-h) / 2.0;
						wert[6] += (-1.0) * 2.0 * h;

						//zeile[vtx_global] += (-1.0)*3.0/2.0*h;
						//zeile[vtx_global-nx] += (-1.0)*(-h)/2.0;
						//zeile[vtx_global-2*nx] += (-1.0)*2.0*h;
					}
				}

				else if (i % nx == 1) //linke Seite, aber nicht vordere/hintere Seite oder Boden/Deckel
				{
					index[0] = vtx_global;
					index[1] = vtx_global + nx;
					index[2] = vtx_global - nx;
					index[3] = vtx_global + nx*ny;
					index[4] = vtx_global - nx*ny;
					index[5] = vtx_global + 1;
					index[6] = vtx_global + 2;

					//zentraler Differenzenquotient ist nur in x-Richtung nicht möglich, deshalb 
					//zuerst normal in y/z-Richtung
					wert[0] = -4.0;
					wert[1] = 1.0;
					wert[2] = 1.0;
					wert[3] = 1.0;
					wert[4] = 1.0;
					//modifizierter Differenzenquotient in x-Richtung
					wert[0] += 11.0 / 38.0;
					wert[5] = -28.0 / 38.0;
					wert[6] = 17.0 / 38.0;

					////zentraler Differenzenquotient ist nur in x-Richtung nicht möglich, deshalb 
					//zuerst normal in y/z-Richtung
					//zeile[vtx_global] = -4.0;
					//zeile[vtx_global+nx] = 1.0;
					//zeile[vtx_global-nx] = 1.0;
					//zeile[vtx_global+nx*ny] = 1.0;
					//zeile[vtx_global-nx*ny] = 1.0;
					////modifizierter Differenzenquotient in x-Richtung
					//zeile[vtx_global] += 11.0/38.0;
					//zeile[vtx_global+1] = -28.0/38.0;
					//zeile[vtx_global+2] = 17.0/38.0;

					//NeumannRB, Normalenvektor ist (1,0,0)
					//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten

					wert[0] += 3.0 / 2.0 * h;
					wert[5] += (-h) / 2.0;
					wert[6] += 2.0 * h;

					//zeile[vtx_global] += 3.0/2.0*h;
					//zeile[vtx_global+1] += (-h)/2.0;
					//zeile[vtx_global+2] += 2.0*h;
				}

				else if (i % nx == 0) //rechte Seite, aber nicht vordere/hintere Seite oder Boden/
				//Deckel
				{
					index[0] = vtx_global;
					index[1] = vtx_global + nx;
					index[2] = vtx_global - nx;
					index[3] = vtx_global + nx*ny;
					index[4] = vtx_global - nx*ny;
					index[5] = vtx_global - 1;
					index[6] = vtx_global - 2;

					//zentraler Differenzenquotient ist nur in x-Richtung nicht möglich, deshalb 
					//zuerst normal in y/z-Richtung
					wert[0] = -4.0;
					wert[1] = 1.0;
					wert[2] = 1.0;
					wert[3] = 1.0;
					wert[4] = 1.0;
					//modifizierter Differenzenquotient in x-Richtung
					wert[0] += 11.0 / 38.0;
					wert[5] = -28.0 / 38.0;
					wert[6] = 17.0 / 38.0;

					////zentraler Differenzenquotient ist nur in x-Richtung nicht möglich, deshalb 
					//zuerst normal in y/z-Richtung
					//zeile[vtx_global] = -4.0;
					//zeile[vtx_global+nx] = 1.0;
					//zeile[vtx_global-nx] = 1.0;
					//zeile[vtx_global+nx*ny] = 1.0;
					//zeile[vtx_global-nx*ny] = 1.0;
					////modifizierter Differenzenquotient in x-Richtung
					//zeile[vtx_global] += 11.0/38.0;
					//zeile[vtx_global-1] = -28.0/38.0;
					//zeile[vtx_global-2] = 17.0/38.0;

					//NeumannRB, Normalenvektor ist (-1,0,0)
					//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten

					wert[0] += (-1.0) * 3.0 / 2.0 * h;
					wert[5] += (-1.0)*(-h) / 2.0;
					wert[6] += (-1.0) * 2.0 * h;

					//zeile[vtx_global] += (-1.0)*3.0/2.0*h;
					//zeile[vtx_global-1] += (-1.0)*(-h)/2.0;
					//zeile[vtx_global-2] += (-1.0)*2.0*h;
				}

				else //innere Punkte
				{
					index[0] = vtx_global;
					index[1] = vtx_global + 1;
					index[2] = vtx_global - 1;
					index[3] = vtx_global + nx;
					index[4] = vtx_global - nx;
					index[5] = vtx_global + nx*ny;
					index[6] = vtx_global - nx*ny;

					//zentraler Differenzenquotient in alle Richtung möglich
					wert[0] = -6.0;
					wert[1] = 1.0;
					wert[2] = 1.0;
					wert[3] = 1.0;
					wert[4] = 1.0;
					wert[5] = 1.0;
					wert[6] = 1.0;

					////zentraler Differenzenquotient in alle Richtung möglich
					//zeile[vtx_global] = -6.0;
					//zeile[vtx_global+1] = 1.0;
					//zeile[vtx_global-1] = 1.0;
					//zeile[vtx_global+nx] = 1.0;
					//zeile[vtx_global-nx] = 1.0;
					//zeile[vtx_global+nx*ny] = 1.0;
					//zeile[vtx_global-nx*ny] = 1.0;

					//keine RB
				}

				if (i!=1) {for (int j=0; j<7; j++) A.sequential_fill(index[j], wert[j]);}
				A.end_of_row();

				rhs.set_local(vtx_global, h*h*bdry(vtx_global));
			}
			return { A, rhs };
		}





















void assembleLeftSidePanel(int** indexMatrix, double** valueMatrix,int fron, int lron, int msize,double h,int nx,int ny, int nz)
{
    int start = 0;
    int end = nx*ny*nz-nx;

    while(start < fron)
    {
        start += nx;
    }

    while(end > lron)
    {
        end -= nx;
    }

    for(int i =start; i<=end ;i+=nx)
    {
        int vtx_local = i - fron;
        int vtx_global = i;
        //fuelle wie linke Seite
        indexMatrix[vtx_local][0] = vtx_global;
		indexMatrix[vtx_local][1] = vtx_global + nx;
		indexMatrix[vtx_local][2] = vtx_global - nx;
		indexMatrix[vtx_local][3] = vtx_global + nx*ny;
		indexMatrix[vtx_local][4] = vtx_global - nx*ny;
		indexMatrix[vtx_local][5] = vtx_global + 1;
		indexMatrix[vtx_local][6] = vtx_global + 2;
		/*
		//zentraler Differenzenquotient ist nur in x-Richtung nicht möglich, deshalb zuerst normal in 
		//y/z-Richtung
		valueMatrix[vtx_local][0] = -4.0;
		valueMatrix[vtx_local][1] = 1.0;
		valueMatrix[vtx_local][2] = 1.0;
		valueMatrix[vtx_local][3] = 1.0;
		valueMatrix[vtx_local][4] = 1.0;
		//modifizierter Differenzenquotient in x-Richtung
		valueMatrix[vtx_local][0] += 1.0;
		valueMatrix[vtx_local][5] = -2.0;
		valueMatrix[vtx_local][6] = 1.0;
		*/
		//NeumannRB, Normalenvektor ist (1,0,0)
		//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
		valueMatrix[vtx_local][0] = 3.0 / 2.0;
		valueMatrix[vtx_local][1] = 0.0;
		valueMatrix[vtx_local][2] = 0.0;
		valueMatrix[vtx_local][3] = 0.0;
		valueMatrix[vtx_local][4] = 0.0;
		valueMatrix[vtx_local][5] = (-4.0) / 2.0;
		valueMatrix[vtx_local][6] = 1.0 / 2.0;
    }
}

void assembleRightSidePanel(int** indexMatrix, double** valueMatrix,int fron, int lron, int msize,double h,int nx,int ny, int nz)
{
    int start = nx-1;
    int end = nx*ny*nz-1;

    while(start < fron)
    {
        start += nx;
    }

    while(end > lron)
    {
        end -= nx;
    }

    for(int i=start;i<=lron;i+=nx)
    {
        int vtx_local = i-fron;
        int vtx_global = i;
        //Fuelle wie rechte seite.
        indexMatrix[vtx_local][0] = vtx_global;
		indexMatrix[vtx_local][1] = vtx_global + nx;
		indexMatrix[vtx_local][2] = vtx_global - nx;
		indexMatrix[vtx_local][3] = vtx_global + nx*ny;
		indexMatrix[vtx_local][4] = vtx_global - nx*ny;
		indexMatrix[vtx_local][5] = vtx_global - 1;
		indexMatrix[vtx_local][6] = vtx_global - 2;
		/*
		//zentraler Differenzenquotient ist nur in x-Richtung nicht möglich, deshalb zuerst normal in 
		//y/z-Richtung
		valueMatrix[vtx_local][0] = -4.0;
		valueMatrix[vtx_local][1] = 1.0;
		valueMatrix[vtx_local][2] = 1.0;
		valueMatrix[vtx_local][3] = 1.0;
		valueMatrix[vtx_local][4] = 1.0;
		//modifizierter Differenzenquotient in x-Richtung
		valueMatrix[vtx_local][0] += 1.0;
		valueMatrix[vtx_local][5] = -2.0;
		valueMatrix[vtx_local][6] = 1.0;
		*/
		//NeumannRB, Normalenvektor ist (-1,0,0)
		//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
		valueMatrix[vtx_local][0] = (-1.0) * 3.0 / 2.0;
		valueMatrix[vtx_local][1] = 0.0;
		valueMatrix[vtx_local][2] = 0.0;
		valueMatrix[vtx_local][3] = 0.0;
		valueMatrix[vtx_local][4] = 0.0;
		valueMatrix[vtx_local][5] = (-1.0) * (-4.0) / 2.0;
		valueMatrix[vtx_local][6] = (-1.0) * 1.0 / 2.0;
    }
}

void assembleTopPanel(int** indexMatrix, double** valueMatrix,int fron, int lron, int msize,double h,int nx,int ny, int nz)
{
    int start = nx*ny*(nz-1);
    int end = nx*ny*nz-1;
    
	if(fron > start)
    {
        start = fron;
    }
	
    if(lron < end)
    {
        end = lron;
    }

    for(int i=start;i<=end;i++)
    {
        int vtx_global =i;
        int vtx_local = vtx_global - fron;
        //Fuelle wie Top Panel
        indexMatrix[vtx_local][0] = vtx_global;
		indexMatrix[vtx_local][1] = vtx_global + 1;
		indexMatrix[vtx_local][2] = vtx_global - 1;
		indexMatrix[vtx_local][3] = vtx_global + nx;
		indexMatrix[vtx_local][4] = vtx_global - nx;
		indexMatrix[vtx_local][5] = vtx_global - nx*ny;
		indexMatrix[vtx_local][6] = vtx_global - 2 * nx*ny;
		/*
		//zentraler Differenzenquotient in x/y-Richtung möglich
		valueMatrix[vtx_local][0] = -4.0;
		valueMatrix[vtx_local][1] = 1.0;
		valueMatrix[vtx_local][2] = 1.0;
		valueMatrix[vtx_local][3] = 1.0;
		valueMatrix[vtx_local][4] = 1.0;
		//modifizierter Differenzenquotient in z-Richtung
		valueMatrix[vtx_local][0] += 1.0;
		valueMatrix[vtx_local][5] = -2.0;
		valueMatrix[vtx_local][6] = 1.0;
		*/
		//NeumannRB, Normalenvektor ist (0,0,-1))
		//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
		valueMatrix[vtx_local][0] = (-1.0) * 3.0 / 2.0;
		valueMatrix[vtx_local][1] = 0.0;
		valueMatrix[vtx_local][2] = 0.0;
		valueMatrix[vtx_local][3] = 0.0;
		valueMatrix[vtx_local][4] = 0.0;
		valueMatrix[vtx_local][5] = (-1.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][6] = (-1.0) * 1.0 / 2.0;
    }
}

void assembleBottomPanel(int** indexMatrix, double** valueMatrix,int fron, int lron, int msize,double h,int nx,int ny, int nz)
{
    int start = 0;
    int end = nx*ny-1;
 
	if(fron > start)
    {
        start = fron;
    }
 
	if(lron < end)
    {
        end = lron;
    }

    for(int i=start;i<=end;i++)
    {
        int vtx_global =i;
        int vtx_local = vtx_global - fron;
		
        indexMatrix[vtx_local][0] = vtx_global;
		indexMatrix[vtx_local][1] = vtx_global + 1;
		indexMatrix[vtx_local][2] = vtx_global - 1;
		indexMatrix[vtx_local][3] = vtx_global + nx;
		indexMatrix[vtx_local][4] = vtx_global - nx;
		indexMatrix[vtx_local][5] = vtx_global + nx*ny;
		indexMatrix[vtx_local][6] = vtx_global + 2 * nx*ny;
		/*
		////zentraler Differenzenquotient in x/y-Richtung möglich
		valueMatrix[vtx_local][0] = -4.0;
		valueMatrix[vtx_local][1] = 1.0;
		valueMatrix[vtx_local][2] = 1.0;
		valueMatrix[vtx_local][3] = 1.0;
		valueMatrix[vtx_local][4] = 1.0;
		//modifizierter Differenzenquotient in z-Richtung
		valueMatrix[vtx_local][0] += 1.0;
		valueMatrix[vtx_local][5] = -2.0;
		valueMatrix[vtx_local][6] = 1.0;
		*/

		//NeumannRB, Normalenvektor ist (0,0,1))
		//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
		valueMatrix[vtx_local][0] = 3.0 / 2.0;
		valueMatrix[vtx_local][1] = 0.0;
		valueMatrix[vtx_local][2] = 0.0;
		valueMatrix[vtx_local][3] = 0.0;
		valueMatrix[vtx_local][4] = 0.0;
		valueMatrix[vtx_local][5] = (-4.0) / 2.0;
		valueMatrix[vtx_local][6] = 1.0 / 2.0;
    }
}

void assembleFrontPanel(int** indexMatrix, double** valueMatrix,int fron, int lron, int msize,double h,int nx,int ny, int nz)
{
    //Hier fuer finde welche punkte im front panel liegen und
    //fuelle alle die im bereich von fron..lron liegen.
    int* frontPanelIdxs = new int[nx*nz];
    for(int iz =0; iz < nz;iz++)
    {
        for(int ix =0;ix<nx;ix++)
        {
            int idx = iz*(nx*ny) + ix;
            frontPanelIdxs[iz*nx+ix] = idx;
        }
    }

    int idxStart =0;
    int start =frontPanelIdxs[idxStart];
    while(start < fron)
    {
        idxStart++;
        start = frontPanelIdxs[idxStart];
    }

    int idxEnd = nx*nz-1;
    int end = frontPanelIdxs[idxEnd];
    while(lron < end)
    {
        idxEnd--;
        end = frontPanelIdxs[idxEnd];
    }

    //Alle Benoetigten idx's liegen zwischen idxStart und idxEnd
    for(int idx = idxStart; idx<=idxEnd;idx++)
    {
        int vtx_global = frontPanelIdxs[idx];
        int vtx_local =vtx_global - fron;
        //fuelle wie in der vorderen seite

        indexMatrix[vtx_local][0] = vtx_global;
		indexMatrix[vtx_local][1] = vtx_global + nx*ny;
		indexMatrix[vtx_local][2] = vtx_global - nx*ny;
		indexMatrix[vtx_local][3] = vtx_global + 1;
		indexMatrix[vtx_local][4] = vtx_global - 1;
		indexMatrix[vtx_local][5] = vtx_global + nx;
		indexMatrix[vtx_local][6] = vtx_global + 2 * nx;
		/*
    	//zentraler Differenzenquotient in x/z-Richtung möglich
		valueMatrix[vtx_local][0] = -4.0;
		valueMatrix[vtx_local][1] = 1.0;
		valueMatrix[vtx_local][2] = 1.0;
		valueMatrix[vtx_local][3] = 1.0;
		valueMatrix[vtx_local][4] = 1.0;
		//modifizierter Differenzenquotient in y-Richtung
		valueMatrix[vtx_local][0] += 1.0;
		valueMatrix[vtx_local][5] = -2.0;
		valueMatrix[vtx_local][6] = 1.0;
		*/
		//NeumannRB, Normalenvektor ist (0,1,0)
	    //RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
		valueMatrix[vtx_local][0] = 3.0 / 2.0;
		valueMatrix[vtx_local][1] = 0.0;
		valueMatrix[vtx_local][2] = 0.0;
		valueMatrix[vtx_local][3] = 0.0;
		valueMatrix[vtx_local][4] = 0.0;
		valueMatrix[vtx_local][5] = (-4.0) / 2.0;
		valueMatrix[vtx_local][6] = 1.0 / 2.0;
    }
}

void assembleBackPanel(int** indexMatrix, double** valueMatrix,int fron, int lron, int msize,double h,int nx,int ny, int nz)
{
    //Hier fuer finde welche punkte im back panel liegen und
    //fuelle alle die im bereich von fron..lron liegen.
    int* backPanelIdxs = new int[nx*nz];
    for(int iz =0; iz < nz;iz++)
    {
        for(int ix =0;ix<nx;ix++)
        {
            int idx = iz*(nx*ny) + ix + nx*(ny-1);
            backPanelIdxs[iz*nx+ix] = idx;
        }
    }
	
    int idxStart =0;
    int start =backPanelIdxs[idxStart];
    while(start < fron)
    {
        idxStart++;
        start = backPanelIdxs[idxStart];
    }

    int idxEnd = nx*nz-1;
    int end = backPanelIdxs[idxEnd];
    while(lron < end)
    {
        idxEnd--;
        end = backPanelIdxs[idxEnd];
    }

    //Alle benoetigten idx's liegen zwischen idxStart und idxEnd
    for(int idx = idxStart; idx<=idxEnd;idx++)
    {
        int vtx_global = backPanelIdxs[idx];
        int vtx_local =vtx_global - fron;
        //fuelle wie in der hinteren Seite
        indexMatrix[vtx_local][0] = vtx_global;
		indexMatrix[vtx_local][1] = vtx_global + nx*ny;
		indexMatrix[vtx_local][2] = vtx_global - nx*ny;
		indexMatrix[vtx_local][3] = vtx_global + 1;
		indexMatrix[vtx_local][4] = vtx_global - 1;
		indexMatrix[vtx_local][5] = vtx_global - nx;
		indexMatrix[vtx_local][6] = vtx_global - 2 * nx;
		/*
		//zentraler Differenzenquotient in x/z-Richtung möglich
		valueMatrix[vtx_local][0] = -4.0;
		valueMatrix[vtx_local][1] = 1.0;
		valueMatrix[vtx_local][2] = 1.0;
		valueMatrix[vtx_local][3] = 1.0;
		valueMatrix[vtx_local][4] = 1.0;
		//modifizierter Differenzenquotient in y-Richtung
		valueMatrix[vtx_local][0] += 1.0;
		valueMatrix[vtx_local][5] = -2.0;
	    valueMatrix[vtx_local][6] = 1.0;
		*/		
		//NeumannRB, Normalenvektor ist (0,-1,0)
		//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
		valueMatrix[vtx_local][0] = (-1.0) * 3.0 / 2.0;
		valueMatrix[vtx_local][1] = 0.0;
		valueMatrix[vtx_local][2] = 0.0;
		valueMatrix[vtx_local][3] = 0.0;
		valueMatrix[vtx_local][4] = 0.0;
		valueMatrix[vtx_local][5] = (-1.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][6] = (-1.0) * 1.0 / 2.0;
    }
}

void assembleKanteVorneLinks(int** indexMatrix, double** valueMatrix,int fron, int lron, int msize,double h,int nx,int ny, int nz)
{
    int* kantenidxs = new int[nz];
    for(int i = 0;i < nz;i++)
    {
        kantenidxs[i] = i*(nx*ny);
    }

    int startidx = 0;
    int endidx = nz;
    while(kantenidxs[startidx] < fron)
    {
        startidx++;
    }
    while(kantenidxs[endidx] > lron)
    {
        endidx--;
    }

    for(int idx = startidx; idx <= endidx ;idx++)
    {
        int vtx_global = kantenidxs[idx];
        int vtx_local = vtx_global - fron;

        indexMatrix[vtx_local][0] = vtx_global;
		indexMatrix[vtx_local][1] = vtx_global + nx*ny;
		indexMatrix[vtx_local][2] = vtx_global - nx*ny;
		indexMatrix[vtx_local][3] = vtx_global + 1;
		indexMatrix[vtx_local][4] = vtx_global + 2;
		indexMatrix[vtx_local][5] = vtx_global + nx;
		indexMatrix[vtx_local][6] = vtx_global + 2 * nx;
		/*
		//zentraler Differenzenquotient nur in z-Richtung möglich
		valueMatrix[vtx_local][0] = -2.0;
		valueMatrix[vtx_local][1] = 1.0;
		valueMatrix[vtx_local][2] = 1.0;
		//modifizierter Differenzenquotient in x/y-Richtung
		valueMatrix[vtx_local][0] += 2.0 * 1.0;
		valueMatrix[vtx_local][3] = -2.0;
		valueMatrix[vtx_local][4] = 1.0;
		valueMatrix[vtx_local][5] = -2.0;
		valueMatrix[vtx_local][6] = 1.0;
		*/
	
		//NeumannRB, Normalenvektor ist (1/sqrt(2),1/sqrt(2),0)
		//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
		valueMatrix[vtx_local][0] = 2.0 * 1.0 / sqrt(2.0) * 3.0 / 2.0;
		valueMatrix[vtx_local][1] = 0.0;
		valueMatrix[vtx_local][2] = 0.0;
		valueMatrix[vtx_local][3] = 1.0 / sqrt(2.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][4] = 1.0 / sqrt(2.0) * 1.0 / 2.0;
		valueMatrix[vtx_local][5] = 1.0 / sqrt(2.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][6] = 1.0 / sqrt(2.0) * 1.0 / 2.0;
    }
}

void assembleKanteVorneRechts(int** indexMatrix, double** valueMatrix,int fron, int lron, int msize,double h,int nx,int ny, int nz)
{
	int* kantenidxs = new int[nz];
    for(int i = 0;i < nz;i++)
    {
        kantenidxs[i] = i*(nx*ny)+(nx-1);
    }

    int startidx = 0;
    int endidx = nz;
    while(kantenidxs[startidx] < fron)
    {
        startidx++;
    }
    while(kantenidxs[endidx] > lron)
    {
        endidx--;
    }
	for(int idx = startidx; idx <= endidx ;idx++)
    {
       int vtx_global = kantenidxs[idx];
       int vtx_local = vtx_global - fron;

        indexMatrix[vtx_local][0] = vtx_global;
		indexMatrix[vtx_local][1] = vtx_global + nx*ny;
		indexMatrix[vtx_local][2] = vtx_global - nx*ny;
		indexMatrix[vtx_local][3] = vtx_global - 1;
		indexMatrix[vtx_local][4] = vtx_global - 2;
		indexMatrix[vtx_local][5] = vtx_global + nx;
		indexMatrix[vtx_local][6] = vtx_global + 2 * nx;
		/*
		//zentraler Differenzenquotient nur in z-Richtung möglich
		valueMatrix[vtx_local][0] = -2.0;
		valueMatrix[vtx_local][1] = 1.0;
		valueMatrix[vtx_local][2] = 1.0;
		//modifizierter Differenzenquotient in x/y-Richtung
		valueMatrix[vtx_local][0] += 2.0 * 1.0;
		valueMatrix[vtx_local][3] = -2.0;
		valueMatrix[vtx_local][4] = 1.0;
		valueMatrix[vtx_local][5] = -2.0;
		valueMatrix[vtx_local][6] = 1.0;
		*/
		//NeumannRB, Normalenvektor ist (-1/sqrt(2),1/sqrt(2),0)
		//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
		valueMatrix[vtx_local][0] = 0.0 * 1.0 / sqrt(2.0) * 3.0 / 2.0;
		valueMatrix[vtx_local][1] = 0.0;
		valueMatrix[vtx_local][2] = 0.0;
		valueMatrix[vtx_local][3] = (-1.0) / sqrt(2.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][4] = (-1.0) / sqrt(2.0) * 1.0 / 2.0;
		valueMatrix[vtx_local][5] = 1.0 / sqrt(2.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][6] = 1.0 / sqrt(2.0) * 1.0 / 2.0;
    }
}

void assembleKanteHintenLinks(int** indexMatrix, double** valueMatrix,int fron, int lron, int msize,double h,int nx,int ny, int nz)
{
	int* kantenidxs = new int[nz];
    for(int i = 0;i < nz;i++)
    {
        kantenidxs[i] = i*(nx*ny)+(nx*ny-1)-(nx-1);
    }

    int startidx = 0;
    int endidx = nz;
    while(kantenidxs[startidx] < fron)
    {
        startidx++;
    }
    while(kantenidxs[endidx] > lron)
    {
        endidx--;
    }

	 for(int idx = startidx; idx <= endidx ;idx++)
    {
        int vtx_global = kantenidxs[idx];
        int vtx_local = vtx_global - fron;

        indexMatrix[vtx_local][0] = vtx_global;
		indexMatrix[vtx_local][1] = vtx_global + nx*ny;
		indexMatrix[vtx_local][2] = vtx_global - nx*ny;
		indexMatrix[vtx_local][3] = vtx_global + 1;
		indexMatrix[vtx_local][4] = vtx_global + 2;
		indexMatrix[vtx_local][5] = vtx_global - nx;
		indexMatrix[vtx_local][6] = vtx_global - 2 * nx;
		/*
		//zentraler Differenzenquotient nur in z-Richtung möglich
		valueMatrix[vtx_local][0] = -2.0;
		valueMatrix[vtx_local][1] = 1.0;
		valueMatrix[vtx_local][2] = 1.0;
		//modifizierter Differenzenquotient in x/y-Richtung
		valueMatrix[vtx_local][0] += 2.0 * 1.0;
		valueMatrix[vtx_local][3] = -2.0;
		valueMatrix[vtx_local][4] = 1.0;
		valueMatrix[vtx_local][5] = -2.0;
		valueMatrix[vtx_local][6] = 1.0;
		*/
		//NeumannRB, Normalenvektor ist (1/sqrt(2),-1/sqrt(2),0)
		//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
		valueMatrix[vtx_local][0] = 0.0 * 1.0 / sqrt(2.0) * 3.0 / 2.0;
		valueMatrix[vtx_local][1] = 0.0;
		valueMatrix[vtx_local][2] = 0.0;
		valueMatrix[vtx_local][3] = 1.0 / sqrt(2.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][4] = 1.0 / sqrt(2.0) * 1.0 / 2.0;
		valueMatrix[vtx_local][5] = (-1.0) / sqrt(2.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][6] = (-1.0) / sqrt(2.0) * 1.0 / 2.0;
    }
}

void assembleKanteHintenRechts(int** indexMatrix, double** valueMatrix,int fron, int lron, int msize,double h,int nx,int ny, int nz)
{
	int* kantenidxs = new int[nz];
    for(int i = 0;i < nz;i++)
    {
        kantenidxs[i] = i*(nx*ny)+(nx*ny-1);
    }

    int startidx = 0;
    int endidx = nz;
    while(kantenidxs[startidx] < fron)
    {
        startidx++;
    }
    while(kantenidxs[endidx] > lron)
    {
        endidx--;
    }

	 for(int idx = startidx; idx <= endidx ;idx++)
    {
        int vtx_global = kantenidxs[idx];
        int vtx_local = vtx_global - fron;

        indexMatrix[vtx_local][0] = vtx_global;
		indexMatrix[vtx_local][1] = vtx_global + nx*ny;
		indexMatrix[vtx_local][2] = vtx_global - nx*ny;
		indexMatrix[vtx_local][3] = vtx_global - 1;
		indexMatrix[vtx_local][4] = vtx_global - 2;
		indexMatrix[vtx_local][5] = vtx_global - nx;
		indexMatrix[vtx_local][6] = vtx_global - 2 * nx;
		/*
		//zentraler Differenzenquotient nur in z-Richtung möglich
		valueMatrix[vtx_local][0] = -2.0;
		valueMatrix[vtx_local][1] = 1.0;
		valueMatrix[vtx_local][2] = 1.0;
		//modifizierter Differenzenquotient in x/y-Richtung
		valueMatrix[vtx_local][0] += 2.0 * 1.0;
		valueMatrix[vtx_local][3] = -2.0;
		valueMatrix[vtx_local][4] = 1.0;
		valueMatrix[vtx_local][5] = -2.0;
		valueMatrix[vtx_local][6] = 1.0;
		*/
		//NeumannRB, Normalenvektor ist (- 1/sqrt(2),- 1/sqrt(2),0)
		//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
		valueMatrix[vtx_local][0] = 2.0 * (-1.0) / sqrt(2.0) * 3.0 / 2.0;
		valueMatrix[vtx_local][1] = 0.0;
		valueMatrix[vtx_local][2] = 0.0;
		valueMatrix[vtx_local][3] = (-1.0) / sqrt(2.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][4] = (-1.0) / sqrt(2.0) * 1.0 / 2.0;
		valueMatrix[vtx_local][5] = (-1.0) / sqrt(2.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][6] = (-1.0) / sqrt(2.0) * 1.0 / 2.0;
    }
}

//Untere Kanten

void assembleKanteUntenVorne(int** indexMatrix, double** valueMatrix,int fron, int lron, int msize,double h,int nx,int ny, int nz)
{
    int* kantenidxs = new int[nx];
    for(int i = 0;i < nx;i++)
    {
        kantenidxs[i] = i;
    }

    int startidx = 0;
    int endidx = nx;
    while(kantenidxs[startidx] < fron)
    {
        startidx++;
    }
    while(kantenidxs[endidx] > lron)
    {
        endidx--;
    }

    for(int idx = startidx; idx <= endidx ;idx++)
    {
        int vtx_global = kantenidxs[idx];
        int vtx_local = vtx_global - fron;
        //Fuelle wie in Kante unten Vorne
        indexMatrix[vtx_local][0] = vtx_global;
		indexMatrix[vtx_local][1] = vtx_global + 1;
		indexMatrix[vtx_local][2] = vtx_global - 1;
		indexMatrix[vtx_local][3] = vtx_global + nx;
		indexMatrix[vtx_local][4] = vtx_global + 2 * nx;
		indexMatrix[vtx_local][5] = vtx_global + nx*ny;
		indexMatrix[vtx_local][6] = vtx_global + 2 * nx*ny;
		/*
		//zentraler Differezenquotient in x-Richtung möglich
		valueMatrix[vtx_local][0] = -2.0;
		valueMatrix[vtx_local][1] = 1.0;
		valueMatrix[vtx_local][2] = 1.0;
		//modifizierter Differenzenquotient in y/z-Richtung
		valueMatrix[vtx_local][0] += 2.0 * 1.0;
		valueMatrix[vtx_local][3] = -2.0;
		valueMatrix[vtx_local][4] = 1.0;
		valueMatrix[vtx_local][5] = -2.0;
		valueMatrix[vtx_local][6] = 1.0;
		*/
		//NeumannRB, Normalenvektor ist (0,1/sqrt(2),1/sqrt(2))
		//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
		valueMatrix[vtx_local][0] = 2.0 * 1.0 / sqrt(2.0) * 3.0 / 2.0;
		valueMatrix[vtx_local][1] = 0.0;
		valueMatrix[vtx_local][2] = 0.0;
		valueMatrix[vtx_local][3] = 1.0 / sqrt(2.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][4] = 1.0 / sqrt(2.0) * 1.0 / 2.0;
		valueMatrix[vtx_local][5] = 1.0 / sqrt(2.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][6] = 1.0 / sqrt(2.0) * 1.0 / 2.0;
    }
}

void assembleKanteUntenHinten(int** indexMatrix, double** valueMatrix,int fron, int lron, int msize,double h,int nx,int ny, int nz)
{
    int* kantenidxs = new int[nx];
    for(int i = 0;i < nx;i++)
    {
        kantenidxs[i] = i+(nx*ny-1)-(nx-1);
    }

    int startidx = 0;
    int endidx = nx;
    while(kantenidxs[startidx] < fron)
    {
        startidx++;
    }
    while(kantenidxs[endidx] > lron)
    {
        endidx--;
    }

    for(int idx = startidx; idx <= endidx ;idx++)
    {
        int vtx_global = kantenidxs[idx];
        int vtx_local = vtx_global - fron;

        indexMatrix[vtx_local][0] = vtx_global;
		indexMatrix[vtx_local][1] = vtx_global + 1;
		indexMatrix[vtx_local][2] = vtx_global - 1;
		indexMatrix[vtx_local][3] = vtx_global - nx;
		indexMatrix[vtx_local][4] = vtx_global - 2 * nx;
		indexMatrix[vtx_local][5] = vtx_global + nx*ny;
		indexMatrix[vtx_local][6] = vtx_global + 2 * nx*ny;
		/*
		//zentraler Differezenquotient in x-Richtung möglich
		valueMatrix[vtx_local][0] = -2.0;
		valueMatrix[vtx_local][1] = 1.0;
		valueMatrix[vtx_local][2] = 1.0;
		//modifizierter Differenzenquotient in y/z-Richtung
		valueMatrix[vtx_local][0] += 2.0 * 1.0;
		valueMatrix[vtx_local][3] = -2.0;
		valueMatrix[vtx_local][4] = 1.0;
		valueMatrix[vtx_local][5] = -2.0;
		valueMatrix[vtx_local][6] = 1.0;
		*/
		
		//NeumannRB, Normalenvektor ist (0,-1/sqrt(2),1/sqrt(2))
		//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
		valueMatrix[vtx_local][0] = 0.0*1.0/sqrt(2.0)*3.0/2.0;
		valueMatrix[vtx_local][1] = 0.0;
		valueMatrix[vtx_local][2] = 0.0;
		valueMatrix[vtx_local][3] = (-1.0) * 1.0 / sqrt(2.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][4] = (-1.0) * 1.0 / sqrt(2.0) * 1.0 / 2.0;
		valueMatrix[vtx_local][5] = 1.0 / sqrt(2.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][6] = 1.0 / sqrt(2.0) * 1.0 / 2.0;
    }
}

void assembleKanteUntenLinks(int** indexMatrix, double** valueMatrix,int fron, int lron, int msize,double h,int nx,int ny, int nz)
{
    int* kantenidxs = new int[ny];
    for(int i = 0;i < ny;i++)
    {
        kantenidxs[i] = i*nx;
    }

    int startidx = 0;
    int endidx = ny;
    while(kantenidxs[startidx] < fron)
    {
        startidx++;
    }
    while(kantenidxs[endidx] > lron)
    {
        endidx--;
    }

    for(int idx = startidx; idx <= endidx ;idx++)
    {
        int vtx_global = kantenidxs[idx];
        int vtx_local = vtx_global - fron;
        indexMatrix[vtx_local][0] = vtx_global;
		indexMatrix[vtx_local][1] = vtx_global + nx;
		indexMatrix[vtx_local][2] = vtx_global - nx;
		indexMatrix[vtx_local][3] = vtx_global + 1;
		indexMatrix[vtx_local][4] = vtx_global + 2;
		indexMatrix[vtx_local][5] = vtx_global + nx*ny;
		indexMatrix[vtx_local][6] = vtx_global + 2 * nx*ny;
		/*
		//zentraler Differezenquotient in y-Richtung möglich
		valueMatrix[vtx_local][0] = -2.0;
		valueMatrix[vtx_local][1] = 1.0;
		valueMatrix[vtx_local][2] = 1.0;
		//modifizierter Differenzenquotient in x/z-Richtung
		valueMatrix[vtx_local][0] += 2.0 * 1.0;
		valueMatrix[vtx_local][3] = -2.0;
		valueMatrix[vtx_local][4] = 1.0;
		valueMatrix[vtx_local][5] = -2.0;
		valueMatrix[vtx_local][6] = 1.0;
		*/

		//NeumannRB, Normalenvektor ist (1/sqrt(2),0,1/sqrt(2))
		//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
		valueMatrix[vtx_local][0] = 2.0 * 1.0 / sqrt(2.0) * 3.0 / 2.0;
		valueMatrix[vtx_local][1] = 0.0;
		valueMatrix[vtx_local][2] = 0.0;
		valueMatrix[vtx_local][3] = 1.0 / sqrt(2.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][4] = 1.0 / sqrt(2.0) * 1.0 / 2.0;
		valueMatrix[vtx_local][5] = 1.0 / sqrt(2.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][6] = 1.0 / sqrt(2.0) * 1.0 / 2.0;

    }
}

void assembleKanteUntenRechts(int** indexMatrix, double** valueMatrix,int fron, int lron, int msize,double h,int nx,int ny, int nz)
{
    int* kantenidxs = new int[ny];
    for(int i = 0;i < ny;i++)
    {
        kantenidxs[i] = i*nx+nx-1;
    }

    int startidx = 0;
    int endidx = ny;
    while(kantenidxs[startidx] < fron)
    {
        startidx++;
    }
    while(kantenidxs[endidx] > lron)
    {
        endidx--;
    }

    for(int idx = startidx; idx <= endidx ;idx++)
    {
        int vtx_global = kantenidxs[idx];
        int vtx_local = vtx_global - fron;
        indexMatrix[vtx_local][0] = vtx_global;
		indexMatrix[vtx_local][1] = vtx_global + nx;
		indexMatrix[vtx_local][2] = vtx_global - nx;
		indexMatrix[vtx_local][3] = vtx_global - 1;
		indexMatrix[vtx_local][4] = vtx_global - 2;
		indexMatrix[vtx_local][5] = vtx_global + nx*ny;
		indexMatrix[vtx_local][6] = vtx_global + 2 * nx*ny;
		/*
		//zentraler Differezenquotient in y-Richtung möglich
		valueMatrix[vtx_local][0] = -2.0;
		valueMatrix[vtx_local][1] = 1.0;
		valueMatrix[vtx_local][2] = 1.0;
		//modifizierter Differenzenquotient in x/z-Richtung
		valueMatrix[vtx_local][0] += 2.0 * 1.0;
		valueMatrix[vtx_local][3] = -2.0;
		valueMatrix[vtx_local][4] = 1.0;
		valueMatrix[vtx_local][5] = -2.0;
		valueMatrix[vtx_local][6] = 1.0;
		*/
		//NeumannRB, Normalenvektor ist (-1/sqrt(2),0,1/sqrt(2))
		//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
		valueMatrix[vtx_local][0] = 0.0*1.0/sqrt(2.0)*3.0/2.0;
		valueMatrix[vtx_local][1] = 0.0;
		valueMatrix[vtx_local][2] = 0.0;
		valueMatrix[vtx_local][3] = (-1.0) * 1.0 / sqrt(2.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][4] = (-1.0) * 1.0 / sqrt(2.0) * 1.0 / 2.0;
		valueMatrix[vtx_local][5] = 1.0 / sqrt(2.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][6] = 1.0 / sqrt(2.0) * 1.0 / 2.0;
    }
}

//Obere Kanten

void assembleKanteObenVorne(int** indexMatrix, double** valueMatrix,int fron, int lron, int msize,double h,int nx,int ny, int nz)
{
    int* kantenidxs = new int[nx];
    for(int i = 0;i < nx;i++)
    {
        kantenidxs[i] = i+(nx*ny*(nz-1));
    }

    int startidx = 0;
    int endidx = nx;
    while(kantenidxs[startidx] < fron)
    {
        startidx++;
    }
    while(kantenidxs[endidx] > lron)
    {
        endidx--;
    }

    for(int idx = startidx; idx <= endidx ;idx++)
    {
        int vtx_global = kantenidxs[idx];
        int vtx_local = vtx_global - fron;
        //Fuelle wie in Kante unten Vorne
        indexMatrix[vtx_local][0] = vtx_global;
		indexMatrix[vtx_local][1] = vtx_global + 1;
		indexMatrix[vtx_local][2] = vtx_global - 1;
		indexMatrix[vtx_local][3] = vtx_global + nx;
		indexMatrix[vtx_local][4] = vtx_global + 2 * nx;
		indexMatrix[vtx_local][5] = vtx_global - nx*ny;
		indexMatrix[vtx_local][6] = vtx_global - 2 * nx*ny;
		/*
		//zentraler Differezenquotient in x-Richtung möglich
		valueMatrix[vtx_local][0] = -2.0;
		valueMatrix[vtx_local][1] = 1.0;
		valueMatrix[vtx_local][2] = 1.0;
		//modifizierter Differenzenquotient in y/z-Richtung
		valueMatrix[vtx_local][0] += 2.0 * 1.0;
		valueMatrix[vtx_local][3] = -2.0;
		valueMatrix[vtx_local][4] = 1.0;
		valueMatrix[vtx_local][5] = -2.0;
		valueMatrix[vtx_local][6] = 1.0;
		*/

		//NeumannRB, Normalenvektor ist (0,1/sqrt(2),-1/sqrt(2))
		//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
		valueMatrix[vtx_local][0] = 0.0 * 1.0 / sqrt(2.0) * 3.0 / 2.0;
		valueMatrix[vtx_local][1] = 0.0;
		valueMatrix[vtx_local][2] = 0.0;
		valueMatrix[vtx_local][3] = 1.0 / sqrt(2.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][4] = 1.0 / sqrt(2.0) * 1.0 / 2.0;
		valueMatrix[vtx_local][5] = -1.0 / sqrt(2.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][6] = -1.0 / sqrt(2.0) * 1.0 / 2.0;
    }
}

void assembleKanteObenHinten(int** indexMatrix, double** valueMatrix,int fron, int lron, int msize,double h,int nx,int ny, int nz)
{
    int* kantenidxs = new int[nx];
    for(int i = 0;i < nx;i++)
    {
        kantenidxs[i] = i+(nx*(ny-1))+(nx*ny*(nz-1));
    }

    int startidx = 0;
    int endidx = nx;
    while(kantenidxs[startidx] < fron)
    {
        startidx++;
    }
    while(kantenidxs[endidx] > lron)
    {
        endidx--;
    }

    for(int idx = startidx; idx <= endidx ;idx++)
    {
        int vtx_global = kantenidxs[idx];
        int vtx_local = vtx_global - fron;

        indexMatrix[vtx_local][0] = vtx_global;
		indexMatrix[vtx_local][1] = vtx_global + 1;
		indexMatrix[vtx_local][2] = vtx_global - 1;
		indexMatrix[vtx_local][3] = vtx_global - nx;
		indexMatrix[vtx_local][4] = vtx_global - 2 * nx;
		indexMatrix[vtx_local][5] = vtx_global - nx*ny;
		indexMatrix[vtx_local][6] = vtx_global - 2 * nx*ny;
		/*
		//zentraler Differezenquotient in x-Richtung möglich
		valueMatrix[vtx_local][0] = -2.0;
		valueMatrix[vtx_local][1] = 1.0;
		valueMatrix[vtx_local][2] = 1.0;
		//modifizierter Differenzenquotient in y/z-Richtung
		valueMatrix[vtx_local][0] += 2.0 * 1.0;
		valueMatrix[vtx_local][3] = -2.0;
		valueMatrix[vtx_local][4] = 1.0;
		valueMatrix[vtx_local][5] = -2.0;
		valueMatrix[vtx_local][6] = 1.0;
		*/
		//NeumannRB, Normalenvektor ist (0,-1/sqrt(2),-1/sqrt(2))
		//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
		valueMatrix[vtx_local][0] = 2.0*(-1.0)/sqrt(2.0)*3.0/2.0;
		valueMatrix[vtx_local][1] = 0.0;
		valueMatrix[vtx_local][2] = 0.0;
		valueMatrix[vtx_local][3] = (-1.0) * 1.0 / sqrt(2.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][4] = (-1.0) * 1.0 / sqrt(2.0) * 1.0 / 2.0;
		valueMatrix[vtx_local][5] = (-1.0) / sqrt(2.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][6] = (-1.0) / sqrt(2.0) * 1.0 / 2.0;
    }
}

void assembleKantenObenLinks(int** indexMatrix, double** valueMatrix,int fron, int lron, int msize,double h,int nx,int ny, int nz)
{
    int* kantenidxs = new int[ny];
    for(int i = 0;i < ny;i++)
    {
        kantenidxs[i] = i*nx+(nx*ny*(nz-1));
    }

    int startidx = 0;
    int endidx = ny;
    while(kantenidxs[startidx] < fron)
    {
        startidx++;
    }
    while(kantenidxs[endidx] > lron)
    {
        endidx--;
    }

    for(int idx = startidx; idx <= endidx ;idx++)
    {
        int vtx_global = kantenidxs[idx];
        int vtx_local = vtx_global - fron;
        indexMatrix[vtx_local][0] = vtx_global;
		indexMatrix[vtx_local][1] = vtx_global + nx;
		indexMatrix[vtx_local][2] = vtx_global - nx;
		indexMatrix[vtx_local][3] = vtx_global + 1;
		indexMatrix[vtx_local][4] = vtx_global + 2;
		indexMatrix[vtx_local][5] = vtx_global - nx*ny;
		indexMatrix[vtx_local][6] = vtx_global - 2 * nx*ny;
		/*
		//zentraler Differezenquotient in y-Richtung möglich
		valueMatrix[vtx_local][0] = -2.0;
		valueMatrix[vtx_local][1] = 1.0;
		valueMatrix[vtx_local][2] = 1.0;
		//modifizierter Differenzenquotient in x/z-Richtung
		valueMatrix[vtx_local][0] += 2.0 * 1.0;
		valueMatrix[vtx_local][3] = -2.0;
		valueMatrix[vtx_local][4] = 1.0;
		valueMatrix[vtx_local][5] = -2.0;
		valueMatrix[vtx_local][6] = 1.0;
		*/

		//NeumannRB, Normalenvektor ist (1/sqrt(2),0,-1/sqrt(2))
		//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
		valueMatrix[vtx_local][0] = 0.0 * 1.0 / sqrt(2.0) * 3.0 / 2.0;
		valueMatrix[vtx_local][1] = 0.0;
		valueMatrix[vtx_local][2] = 0.0;
		valueMatrix[vtx_local][3] = 1.0 / sqrt(2.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][4] = 1.0 / sqrt(2.0) * 1.0 / 2.0;
		valueMatrix[vtx_local][5] = (-1.0) / sqrt(2.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][6] = (-1.0) / sqrt(2.0) * 1.0 / 2.0;
    }
}

void assembleKanteObenRechts(int** indexMatrix, double** valueMatrix,int fron, int lron, int msize,double h,int nx,int ny, int nz)
{
    int* kantenidxs = new int[ny];
    for(int i = 0;i < ny;i++)
    {
        kantenidxs[i] = i*nx+nx-1+(nx*ny*(nz-1));
    }

    int startidx = 0;
    int endidx = ny;
    while(kantenidxs[startidx] < fron)
    {
        startidx++;
    }
    while(kantenidxs[endidx] > lron)
    {
        endidx--;
    }

    for(int idx = startidx; idx <= endidx ;idx++)
    {
        int vtx_global = kantenidxs[idx];
        int vtx_local = vtx_global - fron;
        indexMatrix[vtx_local][0] = vtx_global;
		indexMatrix[vtx_local][1] = vtx_global + nx;
		indexMatrix[vtx_local][2] = vtx_global - nx;
		indexMatrix[vtx_local][3] = vtx_global - 1;
		indexMatrix[vtx_local][4] = vtx_global - 2;
		indexMatrix[vtx_local][5] = vtx_global - nx*ny;
		indexMatrix[vtx_local][6] = vtx_global - 2 * nx*ny;
		/*
		//zentraler Differezenquotient in y-Richtung möglich
		valueMatrix[vtx_local][0] = -2.0;
		valueMatrix[vtx_local][1] = 1.0;
		valueMatrix[vtx_local][2] = 1.0;
		//modifizierter Differenzenquotient in x/z-Richtung
		valueMatrix[vtx_local][0] += 2.0 * 1.0;
		valueMatrix[vtx_local][3] = -2.0;
		valueMatrix[vtx_local][4] = 1.0;
		valueMatrix[vtx_local][5] = -2.0;
		valueMatrix[vtx_local][6] = 1.0;
		*/
	
		//NeumannRB, Normalenvektor ist (-1/sqrt(2),0,-1/sqrt(2))
		//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
		valueMatrix[vtx_local][0] = 2.0*(-1.0)/sqrt(2.0)*3.0/2.0;
		valueMatrix[vtx_local][1] = 0.0;
		valueMatrix[vtx_local][2] = 0.0;
		valueMatrix[vtx_local][3] = (-1.0) * 1.0 / sqrt(2.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][4] = (-1.0) * 1.0 / sqrt(2.0) * 1.0 / 2.0;
		valueMatrix[vtx_local][5] = (-1.0) / sqrt(2.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][6] = (-1.0) / sqrt(2.0) * 1.0 / 2.0;
    }
}

//Die Ecken:
//Die unteren Ecken:
void assembleEckeuntenvornelinks(int** indexMatrix, double** valueMatrix,int fron, int lron, int msize,double h,int nx,int ny, int nz)
{
	//int Eckenidx = 0;
	if(fron <= 0 && 0 <= lron)
	{
		int vtx_global = 0;
		int vtx_local = vtx_global - fron;
	
		indexMatrix[vtx_local][0] = vtx_global;
		indexMatrix[vtx_local][1] = vtx_global + 1;
		indexMatrix[vtx_local][2] = vtx_global + 2;
		indexMatrix[vtx_local][3] = vtx_global + nx;
		indexMatrix[vtx_local][4] = vtx_global + 2*nx;
		indexMatrix[vtx_local][5] = vtx_global + nx*nz;
		indexMatrix[vtx_local][6] = vtx_global + 2 * nx*ny;
	/*
		// Dirichlet-RB
		valueMatrix[vtx_local][0] = 1.0;
		valueMatrix[vtx_local][1] = 0.0;
		valueMatrix[vtx_local][2] = 0.0;
		valueMatrix[vtx_local][3] = 0.0;
		valueMatrix[vtx_local][4] = 0.0;
		valueMatrix[vtx_local][5] = 0.0;
		valueMatrix[vtx_local][6] = 0.0;
		
		//zentraler Differenzenquotient gar nicht möglich
		//modifizierter Differenzenquotient in x/y/z-Richtung
		valueMatrix[vtx_local][0] = 3.0 * 1.0;
		valueMatrix[vtx_local][1] = -2.0;
		valueMatrix[vtx_local][2] = 1.0;
		valueMatrix[vtx_local][3] = -2.0;
		valueMatrix[vtx_local][4] = 1.0;
		valueMatrix[vtx_local][5] = -2.0;
		valueMatrix[vtx_local][6] = 1.0;
	*/
		//NeumannRB, Normalenvektor ist (1/sqrt(3),1/sqrt(3),1/sqrt(3))
		//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
		valueMatrix[vtx_local][0] = 3.0 * 1.0 / sqrt(3.0) * 3.0 / 2.0;
		valueMatrix[vtx_local][1] = 1.0 / sqrt(3.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][2] = 1.0 / sqrt(3.0) * 1.0 / 2.0;
		valueMatrix[vtx_local][3] = 1.0 / sqrt(3.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][4] = 1.0 / sqrt(3.0) * 1.0 / 2.0;
		valueMatrix[vtx_local][5] = 1.0 / sqrt(3.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][6] = 1.0 / sqrt(3.0) * 1.0 / 2.0;
	}
}

void assembleEckeuntenvornerechts(int** indexMatrix, double** valueMatrix,int fron, int lron, int msize,double h,int nx,int ny, int nz)
{
	if(fron <= nx-1 && nx-1 <= lron)
	{
		int vtx_global = nx-1;
		int vtx_local = vtx_global - fron;

		indexMatrix[vtx_local][0] = vtx_global;
		indexMatrix[vtx_local][1] = vtx_global - 1;
		indexMatrix[vtx_local][2] = vtx_global - 2;
		indexMatrix[vtx_local][3] = vtx_global + nx;
		indexMatrix[vtx_local][4] = vtx_global + 2*nx;
		indexMatrix[vtx_local][5] = vtx_global + nx*nz;
		indexMatrix[vtx_local][6] = vtx_global + 2 * nx*ny;
		/*
		//zentraler Differenzenquotient gar nicht möglich
		//modifizierter Differenzenquotient in x/y/z-Richtung
		valueMatrix[vtx_local][0] = 3.0 * 1.0;
		valueMatrix[vtx_local][1] = -2.0;
		valueMatrix[vtx_local][2] = 1.0;
		valueMatrix[vtx_local][3] = -2.0;
		valueMatrix[vtx_local][4] = 1.0;
		valueMatrix[vtx_local][5] = -2.0;
		valueMatrix[vtx_local][6] = 1.0;
		*/
		
		//NeumannRB, Normalenvektor ist (-1/sqrt(3),1/sqrt(3),1/sqrt(3))
		//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
		valueMatrix[vtx_local][0] = 1.0 * 1.0 / sqrt(3.0) * 3.0 / 2.0;
		valueMatrix[vtx_local][1] = (-1.0) / sqrt(3.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][2] = (-1.0) / sqrt(3.0) * 1.0 / 2.0;
		valueMatrix[vtx_local][3] = (1.0) / sqrt(3.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][4] = (1.0) / sqrt(3.0) * 1.0 / 2.0;
		valueMatrix[vtx_local][5] = (1.0) / sqrt(3.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][6] = (1.0) / sqrt(3.0) * 1.0 / 2.0;
	}
}

void assembleEckeuntenhintenlinks(int** indexMatrix, double** valueMatrix,int fron, int lron, int msize,double h,int nx,int ny, int nz)
{
	if(fron <= nx*(ny - 1) && nx*(ny - 1)  <= lron)
	{
		int vtx_global = nx*(ny - 1) ;
		int vtx_local = vtx_global - fron;

		indexMatrix[vtx_local][0] = vtx_global;
		indexMatrix[vtx_local][1] = vtx_global + 1;
		indexMatrix[vtx_local][2] = vtx_global + 2;
		indexMatrix[vtx_local][3] = vtx_global - nx;
		indexMatrix[vtx_local][4] = vtx_global - 2*nx;
		indexMatrix[vtx_local][5] = vtx_global + nx*nz;
		indexMatrix[vtx_local][6] = vtx_global + 2 * nx*ny;
		/*
		//zentraler Differenzenquotient gar nicht möglich
		//modifizierter Differenzenquotient in x/y/z-Richtung
		valueMatrix[vtx_local][0] = 3.0 * 1.0;
		valueMatrix[vtx_local][1] = -2.0;
		valueMatrix[vtx_local][2] = 1.0;
		valueMatrix[vtx_local][3] = -2.0;
		valueMatrix[vtx_local][4] = 1.0;
		valueMatrix[vtx_local][5] = -2.0;
		valueMatrix[vtx_local][6] = 1.0;
		*/

		//NeumannRB, Normalenvektor ist (1/sqrt(3),-1/sqrt(3),1/sqrt(3))
		//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
		valueMatrix[vtx_local][0] = 1.0 * 1.0 / sqrt(3.0) * 3.0 / 2.0;
		valueMatrix[vtx_local][1] = (1.0) / sqrt(3.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][2] = (1.0) / sqrt(3.0) * 1.0 / 2.0;
		valueMatrix[vtx_local][3] = (-1.0) / sqrt(3.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][4] = (-1.0) / sqrt(3.0) * 1.0 / 2.0;
		valueMatrix[vtx_local][5] = (1.0) / sqrt(3.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][6] = (1.0) / sqrt(3.0) * 1.0 / 2.0;
	}
}

void assembleEckeuntenhintenrechts(int** indexMatrix, double** valueMatrix,int fron, int lron, int msize,double h,int nx,int ny, int nz)
{
	if(fron <= nx*ny - 1 && nx*ny - 1 <= lron)
	{
		int vtx_global = nx*ny - 1 ;
		int vtx_local = vtx_global - fron;

		indexMatrix[vtx_local][0] = vtx_global;
		indexMatrix[vtx_local][1] = vtx_global - 1;
		indexMatrix[vtx_local][2] = vtx_global - 2;
		indexMatrix[vtx_local][3] = vtx_global - nx;
		indexMatrix[vtx_local][4] = vtx_global - 2*nx;
		indexMatrix[vtx_local][5] = vtx_global + nx*nz;
		indexMatrix[vtx_local][6] = vtx_global + 2 * nx*ny;
		/*
		//zentraler Differenzenquotient gar nicht möglich
		//modifizierter Differenzenquotient in x/y/z-Richtung
		valueMatrix[vtx_local][0] = 3.0 * 1.0;
		valueMatrix[vtx_local][1] = -2.0;
		valueMatrix[vtx_local][2] = 1.0;
		valueMatrix[vtx_local][3] = -2.0;
		valueMatrix[vtx_local][4] = 1.0;
		valueMatrix[vtx_local][5] = -2.0;
		valueMatrix[vtx_local][6] = 1.0;
		*/
		//NeumannRB, Normalenvektor ist (-1/sqrt(3),-1/sqrt(2),1/sqrt(3))
		//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
		valueMatrix[vtx_local][0] = (-1.0) * 1.0 / sqrt(3.0) * 3.0 / 2.0;
		valueMatrix[vtx_local][1] = (-1.0) / sqrt(3.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][2] = (-1.0) / sqrt(3.0) * 1.0 / 2.0;
		valueMatrix[vtx_local][3] = (-1.0) / sqrt(3.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][4] = (-1.0) / sqrt(3.0) * 1.0 / 2.0;
		valueMatrix[vtx_local][5] = (1.0) / sqrt(3.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][6] = (1.0) / sqrt(3.0) * 1.0 / 2.0;
	}
}

//die oberen ecken:
void assembleEckeObenvornelinks(int** indexMatrix, double** valueMatrix,int fron, int lron, int msize,double h,int nx,int ny, int nz)
{
	int Eckenidx = nx*ny*(nz-1);

	if(fron <= Eckenidx && Eckenidx <= lron)
	{
		int vtx_global = Eckenidx;
		int vtx_local = vtx_global - fron;

        indexMatrix[vtx_local][0] = vtx_global;
		indexMatrix[vtx_local][1] = vtx_global + 1;
		indexMatrix[vtx_local][2] = vtx_global + 2;
		indexMatrix[vtx_local][3] = vtx_global + nx;
		indexMatrix[vtx_local][4] = vtx_global + 2*nx;
		indexMatrix[vtx_local][5] = vtx_global - nx*nz;
		indexMatrix[vtx_local][6] = vtx_global - 2 * nx*ny;
		/*
		//zentraler Differenzenquotient gar nicht möglich
		//modifizierter Differenzenquotient in x/y/z-Richtung
		valueMatrix[vtx_local][0] = 3.0 * 1.0;
		valueMatrix[vtx_local][1] = -2.0;
		valueMatrix[vtx_local][2] = 1.0;
		valueMatrix[vtx_local][3] = -2.0;
		valueMatrix[vtx_local][4] = 1.0;
		valueMatrix[vtx_local][5] = -2.0;
		valueMatrix[vtx_local][6] = 1.0;
		*/

		//NeumannRB, Normalenvektor ist (1/sqrt(3),1/sqrt(3),-1/sqrt(3))
		//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
		valueMatrix[vtx_local][0] = 1.0 * 1.0 / sqrt(3.0) * 3.0 / 2.0;
		valueMatrix[vtx_local][1] = 1.0 / sqrt(3.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][2] = 1.0 / sqrt(3.0) * 1.0 / 2.0;
		valueMatrix[vtx_local][3] = 1.0 / sqrt(3.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][4] = 1.0 / sqrt(3.0) * 1.0 / 2.0;
		valueMatrix[vtx_local][5] = (-1.0) / sqrt(3.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][6] = (-1.0) / sqrt(3.0) * 1.0 / 2.0;
	}
}

void assembleEckeObenvornerechts(int** indexMatrix, double** valueMatrix,int fron, int lron, int msize,double h,int nx,int ny, int nz)
{
    int Eckenidx = nx-1+(nx*ny*(nz-1));
	if(fron <= Eckenidx && Eckenidx <= lron)
	{
		int vtx_global = Eckenidx;
		int vtx_local = vtx_global - fron;

		indexMatrix[vtx_local][0] = vtx_global;
		indexMatrix[vtx_local][1] = vtx_global - 1;
		indexMatrix[vtx_local][2] = vtx_global - 2;
		indexMatrix[vtx_local][3] = vtx_global + nx;
		indexMatrix[vtx_local][4] = vtx_global + 2*nx;
		indexMatrix[vtx_local][5] = vtx_global - nx*nz;
		indexMatrix[vtx_local][6] = vtx_global - 2 * nx*ny;
		/*
		//zentraler Differenzenquotient gar nicht möglich
		//modifizierter Differenzenquotient in x/y/z-Richtung
		valueMatrix[vtx_local][0] = 3.0 * 1.0;
		valueMatrix[vtx_local][1] = -2.0;
		valueMatrix[vtx_local][2] = 1.0;
		valueMatrix[vtx_local][3] = -2.0;
		valueMatrix[vtx_local][4] = 1.0;
		valueMatrix[vtx_local][5] = -2.0;
		valueMatrix[vtx_local][6] = 1.0;
		*/

		//NeumannRB, Normalenvektor ist (-1/sqrt(3),1/sqrt(3),-1/sqrt(3))
		//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
		valueMatrix[vtx_local][0] = (-1.0) * 1.0 / sqrt(3.0) * 3.0 / 2.0;
		valueMatrix[vtx_local][1] = (-1.0) / sqrt(3.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][2] = (-1.0) / sqrt(3.0) * 1.0 / 2.0;
		valueMatrix[vtx_local][3] = (1.0) / sqrt(3.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][4] = (1.0) / sqrt(3.0) * 1.0 / 2.0;
		valueMatrix[vtx_local][5] = (-1.0) / sqrt(3.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][6] = (-1.0) / sqrt(3.0) * 1.0 / 2.0;
	}
}

void assembleEckeObenhintenlinks(int** indexMatrix, double** valueMatrix,int fron, int lron, int msize,double h,int nx,int ny, int nz)
{
    int Eckenidx = nx*(ny - 1) +nx*ny*(nz-1);
	if(fron <= Eckenidx && Eckenidx  <= lron)
	{
		int vtx_global = Eckenidx;
		int vtx_local = vtx_global - fron;

		indexMatrix[vtx_local][0] = vtx_global;
		indexMatrix[vtx_local][1] = vtx_global + 1;
		indexMatrix[vtx_local][2] = vtx_global + 2;
		indexMatrix[vtx_local][3] = vtx_global - nx;
		indexMatrix[vtx_local][4] = vtx_global - 2*nx;
		indexMatrix[vtx_local][5] = vtx_global - nx*nz;
		indexMatrix[vtx_local][6] = vtx_global - 2 * nx*ny;
		/*
		//zentraler Differenzenquotient gar nicht möglich
		//modifizierter Differenzenquotient in x/y/z-Richtung
		valueMatrix[vtx_local][0] = 3.0 * 1.0;
		valueMatrix[vtx_local][1] = -2.0;
		valueMatrix[vtx_local][2] = 1.0;
		valueMatrix[vtx_local][3] = -2.0;
		valueMatrix[vtx_local][4] = 1.0;
		valueMatrix[vtx_local][5] = -2.0;
		valueMatrix[vtx_local][6] = 1.0;
		*/

		//NeumannRB, Normalenvektor ist (1/sqrt(3),-1/sqrt(3),-1/sqrt(3))
		//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten

		valueMatrix[vtx_local][0] = (-1.0) * 1.0 / sqrt(3.0) * 3.0 / 2.0;
		valueMatrix[vtx_local][1] = 1.0 / sqrt(3.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][2] = 1.0 / sqrt(3.0) * 1.0 / 2.0;
		valueMatrix[vtx_local][3] = (-1.0) / sqrt(3.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][4] = (-1.0) / sqrt(3.0) * 1.0 / 2.0;
		valueMatrix[vtx_local][5] = (-1.0) / sqrt(3.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][6] = (-1.0) / sqrt(3.0) * 1.0 / 2.0;
	}
}

void assembleEckeObenhintenrechts(int** indexMatrix, double** valueMatrix,int fron, int lron, int msize,double h,int nx,int ny, int nz)
{
    int Eckenidx = nx*ny*nz-1;
	if(fron <= Eckenidx && Eckenidx <= lron)
	{
		int vtx_global = Eckenidx;
		int vtx_local = vtx_global - fron;

		indexMatrix[vtx_local][0] = vtx_global;
		indexMatrix[vtx_local][1] = vtx_global - 1;
		indexMatrix[vtx_local][2] = vtx_global - 2;
		indexMatrix[vtx_local][3] = vtx_global - nx;
		indexMatrix[vtx_local][4] = vtx_global - 2*nx;
		indexMatrix[vtx_local][5] = vtx_global - nx*nz;
		indexMatrix[vtx_local][6] = vtx_global - 2 * nx*ny;
		/*
		//zentraler Differenzenquotient gar nicht möglich
		//modifizierter Differenzenquotient in x/y/z-Richtung
		valueMatrix[vtx_local][0] = 3.0 * 1.0;
		valueMatrix[vtx_local][1] = -2.0;
		valueMatrix[vtx_local][2] = 1.0;
		valueMatrix[vtx_local][3] = -2.0;
		valueMatrix[vtx_local][4] = 1.0;
		valueMatrix[vtx_local][5] = -2.0;
		valueMatrix[vtx_local][6] = 1.0;
		*/
		//NeumannRB, Normalenvektor ist (-1/sqrt(3),-1/sqrt(3),-1/sqrt(3))
		//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
		valueMatrix[vtx_local][0] = (3.0) * (-1.0) / sqrt(3.0) * 3.0 / 2.0;
		valueMatrix[vtx_local][1] = (-1.0) / sqrt(3.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][2] = (-1.0) / sqrt(3.0) * 1.0 / 2.0;
		valueMatrix[vtx_local][3] = (-1.0) / sqrt(3.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][4] = (-1.0) / sqrt(3.0) * 1.0 / 2.0;
		valueMatrix[vtx_local][5] = (-1.0) / sqrt(3.0)*(-4.0) / 2.0;
		valueMatrix[vtx_local][6] = (-1.0) / sqrt(3.0) * 1.0 / 2.0;
	}
}

template<typename Scalar>
std::pair < DistEllpackMatrix<Scalar>, SlicedVector < Scalar >>
    assemble_neumann_unrolled(size_t nx, size_t ny, size_t nz,
        typename ScalarTraits<Scalar>::RealType h,
        std::function<Scalar(size_t)> bdry)
{
    const size_t N = nx*ny*nz;
    DistEllpackMatrix<Scalar> A(N);
    SlicedVector<Scalar> rhs(N);

    size_t fron = A.first_row_on_node();
    size_t lron = fron + A.get_dim_local() - 1;



    int** indexMatrix = new int*[(int)A.get_dim_local()];
    double** valueMatrix = new double*[(int)A.get_dim_local()];
    for(int i=0;i< A.get_dim_local();i++)
    {
        indexMatrix[i] = new int[7];
        valueMatrix[i] = new double[7];
    }
    int msize = A.get_dim_local();

    //Neuer Plan: Fuelle die Matrix zunaechst als inneres und ueberschreibe
    //danach die Seiten danach Kanten danach Ecken.
    //Das Innere
    for(int i =0;i<A.get_dim_local();i++)
    {
        int vtx_global = i + fron;
        indexMatrix[i][0] = vtx_global;
        indexMatrix[i][1] = vtx_global + 1;
        indexMatrix[i][2] = vtx_global - 1;
        indexMatrix[i][3] = vtx_global + nx;
        indexMatrix[i][4] = vtx_global - nx;
        indexMatrix[i][5] = vtx_global + nx*ny;
        indexMatrix[i][6] = vtx_global - nx*ny;

        //zentraler Differenzenquotient in alle Richtung moeglich
        valueMatrix[i][0] = -6.0;
        valueMatrix[i][1] = 1.0;
        valueMatrix[i][2] = 1.0;
        valueMatrix[i][3] = 1.0;
        valueMatrix[i][4] = 1.0;
        valueMatrix[i][5] = 1.0;
        valueMatrix[i][6] = 1.0;

        ////zentraler Differenzenquotient in alle Richtung moeglich
        //keine RB
    }

    //Fuelle nun die Seiten:
    //Links
    assembleLeftSidePanel(indexMatrix,valueMatrix,fron,lron,msize,h,nx,ny,nz);
    //Rechts
    assembleRightSidePanel(indexMatrix,valueMatrix,fron,lron,msize,h,nx,ny,nz);
    //oben
    assembleTopPanel(indexMatrix,valueMatrix,fron,lron,msize,h,nx,ny,nz);
    //unten
    assembleBottomPanel(indexMatrix,valueMatrix,fron,lron,msize,h,nx,ny,nz);
    //vorne
    assembleFrontPanel(indexMatrix,valueMatrix,fron,lron,msize,h,nx,ny,nz);
    //hinten
    assembleBackPanel(indexMatrix,valueMatrix,fron,lron,msize,h,nx,ny,nz);

    //Fuelle nun die Kanten:
    //untere Kanten
    assembleKanteUntenHinten(indexMatrix,valueMatrix,fron,lron,msize,h,nx,ny,nz);
    assembleKanteUntenRechts(indexMatrix,valueMatrix,fron,lron,msize,h,nx,ny,nz);
    assembleKanteUntenVorne(indexMatrix,valueMatrix,fron,lron,msize,h,nx,ny,nz);
    assembleKanteUntenLinks(indexMatrix,valueMatrix,fron,lron,msize,h,nx,ny,nz);
    //mittlere Kanten
    assembleKanteHintenLinks(indexMatrix,valueMatrix,fron,lron,msize,h,nx,ny,nz);
    assembleKanteHintenRechts(indexMatrix,valueMatrix,fron,lron,msize,h,nx,ny,nz);
    assembleKanteVorneLinks(indexMatrix,valueMatrix,fron,lron,msize,h,nx,ny,nz);
    assembleKanteVorneRechts(indexMatrix,valueMatrix,fron,lron,msize,h,nx,ny,nz);
    //obere Kanten
    assembleKanteObenHinten(indexMatrix,valueMatrix,fron,lron,msize,h,nx,ny,nz);
    assembleKanteObenRechts(indexMatrix,valueMatrix,fron,lron,msize,h,nx,ny,nz);
    assembleKantenObenLinks(indexMatrix,valueMatrix,fron,lron,msize,h,nx,ny,nz);
    assembleKanteObenVorne(indexMatrix,valueMatrix,fron,lron,msize,h,nx,ny,nz);

    //Fuelle nun die Ecken:
    //untere ecken:
    assembleEckeuntenhintenlinks(indexMatrix,valueMatrix,fron,lron,msize,h,nx,ny,nz);
    assembleEckeuntenhintenrechts(indexMatrix,valueMatrix,fron,lron,msize,h,nx,ny,nz);
    assembleEckeuntenvornelinks(indexMatrix,valueMatrix,fron,lron,msize,h,nx,ny,nz);
    assembleEckeuntenvornerechts(indexMatrix,valueMatrix,fron,lron,msize,h,nx,ny,nz);
    //obere ecken:
    assembleEckeObenvornerechts(indexMatrix,valueMatrix,fron,lron,msize,h,nx,ny,nz);
    assembleEckeObenvornelinks(indexMatrix,valueMatrix,fron,lron,msize,h,nx,ny,nz);
    assembleEckeObenhintenlinks(indexMatrix,valueMatrix,fron,lron,msize,h,nx,ny,nz);
    assembleEckeObenhintenrechts(indexMatrix,valueMatrix,fron,lron,msize,h,nx,ny,nz);

    A.prepare_sequential_fill(7);
    
    for(int i =0;i< msize;i++)
    {
        for(int j = 0; j<7;j++)
        {

            A.sequential_fill(indexMatrix[i][j],valueMatrix[i][j]);
        }
        A.end_of_row();
        int vtx_global = indexMatrix[i][0];

        rhs.set_local(vtx_global, h*bdry(vtx_global));
    }
    return {A,rhs};
}

}
//namespace Icarus
