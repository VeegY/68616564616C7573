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
						//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten

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
						//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten

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
						//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten

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
						//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten

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
						//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten

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
						//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten

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
						//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten

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
						//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten

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
						//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten

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
						//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten

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
						//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten

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
						//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten

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
						//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten

						wert[0] += (-1.0) * 3.0 / 2.0 * h;
						wert[5] += (-1.0)*(-h) / 2.0;
						wert[6] += (-1.0) * 2.0 * h;

						//zeile[vtx_global] += (-1.0)*3.0/2.0*h;
						//zeile[vtx_global-nx*ny] += (-1.0)*(-h)/2.0;
						//zeile[vtx_global-2*nx*ny] += (-1.0)*2.0*h;

					}
				}

				else if (i % nx*ny <= nx) //vordere Seite, aber nicht Boden oder Deckel
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
						//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten

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
						//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten

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
						//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten

						wert[0] += 3.0 / 2.0 * h;
						wert[5] += (-h) / 2.0;
						wert[6] += 2.0 * h;

						//zeile[vtx_global] += 3.0/2.0*h;
						//zeile[vtx_global+nx] += (-h)/2.0;
						//zeile[vtx_global+2*nx] += 2.0*h;
					}
				}

				else if (i % nx*ny > nx*(ny - 1)) //hintere Seite, aber nicht Boden oder Deckel
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
						//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten

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
						//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
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
						//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
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

					//zentraler Differenzenquotient ist nur in x-Richtung nicht möglich, deshalb zuerst normal in y/z-Richtung
					wert[0] = -4.0;
					wert[1] = 1.0;
					wert[2] = 1.0;
					wert[3] = 1.0;
					wert[4] = 1.0;
					//modifizierter Differenzenquotient in x-Richtung
					wert[0] += 11.0 / 38.0;
					wert[5] = -28.0 / 38.0;
					wert[6] = 17.0 / 38.0;

					////zentraler Differenzenquotient ist nur in x-Richtung nicht möglich, deshalb zuerst normal in y/z-Richtung
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

				else if (i % nx == 0) //rechte Seite, aber nicht vordere/hintere Seite oder Boden/Deckel
				{
					index[0] = vtx_global;
					index[1] = vtx_global + nx;
					index[2] = vtx_global - nx;
					index[3] = vtx_global + nx*ny;
					index[4] = vtx_global - nx*ny;
					index[5] = vtx_global - 1;
					index[6] = vtx_global - 2;

					//zentraler Differenzenquotient ist nur in x-Richtung nicht möglich, deshalb zuerst normal in y/z-Richtung
					wert[0] = -4.0;
					wert[1] = 1.0;
					wert[2] = 1.0;
					wert[3] = 1.0;
					wert[4] = 1.0;
					//modifizierter Differenzenquotient in x-Richtung
					wert[0] += 11.0 / 38.0;
					wert[5] = -28.0 / 38.0;
					wert[6] = 17.0 / 38.0;

					////zentraler Differenzenquotient ist nur in x-Richtung nicht möglich, deshalb zuerst normal in y/z-Richtung
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
}
