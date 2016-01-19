#ifndef __ASSEMBLE_HPP_
#define __ASSEMBLE_HPP_

#include <functional>

#include "mpihandler.hpp"
#include "distellpackmatrix.hpp"
#include "slicedvector.hpp"
#include "scalartraits.hpp"
#include "utility.hpp"

namespace Icarus
{
/**
 * \brief	Assembliert eine Matrix	mittels Differenzenquotienten zweiter Ordnung.
 *		Es werden Neumann BC verwendet
 *
 * \param A 			Systemmatrix 
 * \param nx			Anzahl der (äquidistanten Punkte in x-Richtung)
 * \param ny			Anzahl der (äquidistanten Punkte in y-Richtung)
 * \param nz			Anzahl der (äquidistanten Punkte in z-Richtung)
 * \param N			Anzahl aller Gitterpunkte
 * \param h			Schrittweite des finiten Differnezenquotienten
 * \function bdry		Funktion, die den Neumann-Wert eines Punktes zurückgibt.
 *
 * \return 			Gibt eine Matrix zurück.
 */
template<typename Scalar>
std::pair<DistEllpackMatrix<Scalar>
assemble_neumann(size_t nx, size_t ny, size_t nz,
		typename ScalarTraits<Scalar>::RealType h,
		std::function<Scalar(size_t)> bdry)
{
    const size_t N = nx*ny*nz;
    DistEllpackMatrix<Scalar> A(N);

    size_t fron = A.first_row_on_node();
    size_t lron = fron + A.get_dim_local() - 1;

	A.prepare_sequential_fill(7);

    for(size_t vtx_global=fron; vtx_global<=lron; vtx_global++)
    {
        std::vector<int> index = { 0, 0, 0, 0, 0, 0, 0 };
        std::vector<double> wert = { 0, 0, 0, 0, 0, 0, 0 };
		
		//vereinfacht den Zugriff auf die Array-Elemente
		int j=size_t;
	
		//Überprüfung der Position
		if(i % nx*ny*nz <= nx*ny)) //Boden
		{
			if(i == 1) //vorderer unterer linker Eckpunkt
			{				
				index[0]=j;
				index[1]=j+1;
				index[2]=j+2;
				index[3]=j+nx;
				index[4]=j+2*nx;
				index[5]=j+nx*ny;
				index[6]=j+2*nx*ny;
				
				//zentraler Differenzenquotient in keine Richtung möglich
				
				wert[0] = 3*11/38;
				wert[1] = -28/38;
				wert[2] = 17/38;
				wert[3] = -28/38;
				wert[4] = 17/38;
				wert[5] = -28/38;
				wert[6] = 17/38;
				
				//zeile[j] = 3*11/38;
				//zeile[j+1] = -28/38;
				//zeile[j+2] = 17/38;
				//zeile[j+nx] = -28/38;
				//zeile[j+2*nx] = 17/38;
				//zeile[j+nx*ny] = -28/38;
				//zeile[j+2*nx*ny] = 17/38;
			
				//NeumannRB, Normalenvektor ist (1/sqrt(3),1/sqrt(3),1/sqrt(3))
				//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
				
				wert[0] += 3*1/sqrt(3)*3/2*h;
				wert[1] += 1/sqrt(3)*(-h)/2;
				wert[2] += 1/sqrt(3) *2*h;
				wert[3] += 1/sqrt(3)*(-h)/2;
				wert[4] += 1/sqrt(3) *2*h;
				wert[5] += 1/sqrt(3)*(-h)/2;
				wert[6] += 1/sqrt(3) *2*h;
				
				//zeile[j] += 3*1/sqrt(3)*3/2*h;
				//zeile[j+1] += 1/sqrt(3)*(-h)/2;
				//zeile[j+2] += 1/sqrt(3) *2*h;
				//zeile[j+nx] += 1/sqrt(3)*(-h)/2;
				//zeile[j+2*nx] += 1/sqrt(3) *2*h;
				//zeile[j+nx*ny] += 1/sqrt(3)*(-h)/2;
				//zeile[j+2*nx*ny] += 1/sqrt(3) *2*h;
			
			}
		
			else if(i == nx) //vorderer unterer rechter Eckpunkt
			{			
				index[0]=j;
				index[1]=j-1;
				index[2]=j-2;
				index[3]=j+nx;
				index[4]=j+2*nx;
				index[5]=j+nx*ny;
				index[6]=j+2*nx*ny;
				
				//zentraler Differenzenquotient in keine Richtung möglich
				
				wert[0] = 3*11/38;
				wert[1] = -28/38;
				wert[2] = 17/38;
				wert[3] = -28/38;
				wert[4] = 17/38;
				wert[5] = -28/38;
				wert[6] = 17/38;
				
				
				//zeile[j] = 3*11/38;
				//zeile[j-1] = -28/38;
				//zeile[j-2] = 17/38;
				//zeile[j+nx] = -28/38;
				//zeile[j+2*nx] = 17/38;
				//zeile[j+nx*ny] = -28/38;
				//zeile[j+2*nx*ny] = 17/38;
			
				//NeumannRB, Normalenvektor ist (-1/sqrt(3),1/sqrt(3),1/sqrt(3))
				//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
				
				wert[0] += 1*1/sqrt(3)*3/2*h;
				wert[1] += (-1)*1/sqrt(3)*(-h)/2;
				wert[2] += (-1)*1/sqrt(3) *2*h;
				wert[3] += 1/sqrt(3)*(-h)/2;
				wert[4] += 1/sqrt(3) *2*h;
				wert[5] += 1/sqrt(3)*(-h)/2;
				wert[6] += 1/sqrt(3) *2*h;
				
				//zeile[j] += 1*1/sqrt(3)*3/2*h;
				//zeile[j-1] += (-1)*1/sqrt(3)*(-h)/2;
				//zeile[j-2] += (-1)*1/sqrt(3) *2*h;
				//zeile[j+nx] += 1/sqrt(3)*(-h)/2;
				//zeile[j+2*nx] += 1/sqrt(3) *2*h;
				//zeile[j+nx*ny] += 1/sqrt(3)*(-h)/2;
				//zeile[j+2*nx*ny] += 1/sqrt(3) *2*h;
			}
		
			else if(i == nx*(ny-1)+1) //hinterer unterer linker Eckpunkt
			{
				index[0]=j;
				index[1]=j+1;
				index[2]=j+2;
				index[3]=j-nx;
				index[4]=j-2*nx;
				index[5]=j+nx*ny;
				index[6]=j+2*nx*ny;
				
				//zentraler Differenzenquotient in keine Richtung möglich
				
				wert[0] = 3*11/38;
				wert[1] = -28/38;
				wert[2] = 17/38;
				wert[3] = -28/38;
				wert[4] = 17/38;
				wert[5] = -28/38;
				wert[6] = 17/38;
				
				//zeile[j] = 3*11/38;
				//zeile[j+1] = -28/38;
				//zeile[j+2] = 17/38;
				//zeile[j-nx] = -28/38;
				//zeile[j-2*nx] = 17/38;
				//zeile[j+nx*ny] = -28/38;
				//zeile[j+2*nx*ny] = 17/38;
			
				//NeumannRB, Normalenvektor ist (1/sqrt(3),-1/sqrt(3),1/sqrt(3))
				//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
				
				wert[0] += 1*1/sqrt(3)*3/2*h;
				wert[1] += 1/sqrt(3)*(-h)/2;
				wert[2] += 1/sqrt(3) *2*h;
				wert[3] += (-1)*1/sqrt(3)*(-h)/2;
				wert[4] += (-1)*1/sqrt(3) *2*h;
				wert[5] += 1/sqrt(3)*(-h)/2;
				wert[6] += 1/sqrt(3) *2*h;
				
				//zeile[j] += 1*1/sqrt(3)*3/2*h;
				//zeile[j+1] += 1/sqrt(3)*(-h)/2;
				//zeile[j+2] += 1/sqrt(3) *2*h;
				//zeile[j-nx] += (-1)*1/sqrt(3)*(-h)/2;
				//zeile[j-2*nx] += (-1)*1/sqrt(3) *2*h;
				//zeile[j+nx*ny] += 1/sqrt(3)*(-h)/2;
				//zeile[j+2*nx*ny] += 1/sqrt(3) *2*h;
			}
		
			else if(i == nx*ny) //hinterer  unterer rechter Eckpunkt
			{	
				index[0]=j;
				index[1]=j-1;
				index[2]=j-2;
				index[3]=j-nx;
				index[4]=j-2*nx;
				index[5]=j+nx*ny;
				index[6]=j+2*nx*ny;
				
				//zentraler Differenzenquotient in keine Richtung möglich
				
				wert[0] = 3*11/38;
				wert[1] = -28/38;
				wert[2] = 17/38;
				wert[3] = -28/38;
				wert[4] = 17/38;
				wert[5] = -28/38;
				wert[6] = 17/38;
				
				//zeile[j] = 3*11/38;
				//zeile[j-1] = -28/38;
				//zeile[j-2] = 17/38;
				//zeile[j-nx] = -28/38;
				//zeile[j-2*nx] = 17/38;
				//zeile[j+nx*ny] = -28/38;
				//zeile[j+2*nx*ny] = 17/38;
			
				//NeumannRB, Normalenvektor ist (-1/sqrt(3),-1/sqrt(3),1/sqrt(3))
				//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
				
				wert[0] += (-1)*1/sqrt(3)*3/2*h;
				wert[1] += (-1)*1/sqrt(3)*(-h)/2;
				wert[2] += (-1)*1/sqrt(3) *2*h;
				wert[3] += (-1)*1/sqrt(3)*(-h)/2;
				wert[4] += (-1)*1/sqrt(3) *2*h;
				wert[5] += 1/sqrt(3)*(-h)/2;
				wert[6] += 1/sqrt(3) *2*h;
				
				//zeile[j] += (-1)*1/sqrt(3)*3/2*h;
				//zeile[j-1] += (-1)*1/sqrt(3)*(-h)/2;
				//zeile[j-2] += (-1)*1/sqrt(3) *2*h;
				//zeile[j-nx] += (-1)*1/sqrt(3)*(-h)/2;
				//zeile[j-2*nx] += (-1)*1/sqrt(3) *2*h;
				//zeile[j+nx*ny] += 1/sqrt(3)*(-h)/2;
				//zeile[j+2*nx*ny] += 1/sqrt(3) *2*h;
			}
		
			else if(i % nx <= nx) //vordere untere Kante ohne Eckpunkt
			{
				index[0]=j;
				index[1]=j+1;
				index[2]=j-1;
				index[3]=j+nx;
				index[4]=j+2*nx;
				index[5]=j+nx*ny;
				index[6]=j+2*nx*ny;
				
				//zentraler Differezenquotient in x-Richtung möglich
				wert[0] = -2;
				wert[1] = 1;
				wert[2] = 1;
				//modifizierter Differenzenquotient in y/z-Richtung
				wert[0] += 2*11/38;
				wert[3] = -28/38;
				wert[4] = 17/38;
				wert[5] = -28/38;
				wert[6] = 17/38;
				
				////zentraler Differezenquotient in x-Richtung möglich
				//zeile[j] = -2;
				//zeile[j+1] = 1;
				//zeile[j-1] = 1;
				////modifizierter Differenzenquotient in y/z-Richtung
				//zeile[j] += 2*11/38;
				//zeile[j+nx] = -28/38;
				//zeile[j+2*nx] = 17/38;
				//zeile[j+nx*ny] = -28/38;
				//zeile[j+2*nx*ny] = 17/38;
				
				//NeumannRB, Normalenvektor ist (0,1/sqrt(2),1/sqrt(2))
				//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
				
				wert[0] += 2*1/sqrt(3)*3/2*h;
				wert[3] += 1/sqrt(3)*(-h)/2;
				wert[4] += 1/sqrt(3) *2*h;
				wert[5] += 1/sqrt(3)*(-h)/2;
				wert[6] += 1/sqrt(3) *2*h;
				
				//zeile[j] += 2*1/sqrt(3)*3/2*h;
				//zeile[j+nx] += 1/sqrt(3)*(-h)/2;
				//zeile[j+2*nx] += 1/sqrt(3) *2*h;
				//zeile[j+nx*ny] += 1/sqrt(3)*(-h)/2;
				//zeile[j+2*nx*ny] += 1/sqrt(3) *2*h;
			}
		
			else if(i % nx*ny > nx*(ny-1)) //hintere untere Kante ohne Eckpunkt
			{
				index[0]=j;
				index[1]=j+1;
				index[2]=j-1;
				index[3]=j-nx;
				index[4]=j-2*nx;
				index[5]=j+nx*ny;
				index[6]=j+2*nx*ny;
				
				//zentraler Differezenquotient in x-Richtung möglich
				wert[0] = -2;
				wert[1] = 1;
				wert[2] = 1;
				//modifizierter Differenzenquotient in y/z-Richtung
				wert[0] += 2*11/38;
				wert[3] = -28/38;
				wert[4] = 17/38;
				wert[5] = -28/38;
				wert[6] = 17/38;
				
				////zentraler Differezenquotient in x-Richtung möglich
				//zeile[j] = -2;
				//zeile[j+1] = 1;
				//zeile[j-1] = 1;
				//modifizierter Differenzenquotient in y/z-Richtung
				//zeile[j] += 2*11/38;
				//zeile[j-nx] = -28/38;
				//zeile[j-2*nx] = 17/38;
				//zeile[j+nx*ny] = -28/38;
				//zeile[j+2*nx*ny] = 17/38;
			
				//NeumannRB, Normalenvektor ist (0,-1/sqrt(2),1/sqrt(2))
				//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
				
				//wert[0] += 0*1/sqrt(3)*3/2*h;
				wert[3] += (-1)*1/sqrt(3)*(-h)/2;
				wert[4] += (-1)*1/sqrt(3) *2*h;
				wert[5] += 1/sqrt(3)*(-h)/2;
				wert[6] += 1/sqrt(3) *2*h;
				
				//zeile[j] += 0*1/sqrt(3)*3/2*h;		
				//zeile[j-nx] += (-1)*1/sqrt(3)*(-h)/2;
				//zeile[j-2*nx] += (-1)*1/sqrt(3) *2*h;
				//zeile[j+nx*ny] += 1/sqrt(3)*(-h)/2;
				//zeile[j+2*nx*ny] += 1/sqrt(3) *2*h;
			}
		
			else if(i % nx == 1) //linke untere Kante ohne Eckpunkt
			{
				index[0]=j;
				index[1]=j+nx;
				index[2]=j-nx;
				index[3]=j+1;
				index[4]=j+2;
				index[5]=j+nx*ny;
				index[6]=j+2*nx*ny;
				
				//zentraler Differezenquotient in y-Richtung möglich
				wert[0] = -2;
				wert[1] = 1;
				wert[2] = 1;
				//modifizierter Differenzenquotient in x/z-Richtung
				wert[0] += 2*11/38;
				wert[3] = -28/38;
				wert[4] = 17/38;
				wert[5] = -28/38;
				wert[6] = 17/38;
				
				////zentraler Differezenquotient in y-Richtung möglich
				//zeile[j] = -2;
				//zeile[j+nx] = 1;
				//zeile[j-nx] = 1;
				////modifizierter Differenzenquotient in x/z-Richtung
				//zeile[j] += 2*11/38;
				//zeile[j+1] = -28/38;
				//zeile[j+2] = 17/38;
				//zeile[j+nx*ny] = -28/38;
				//zeile[j+2*nx*ny] = 17/38;
				
				//NeumannRB, Normalenvektor ist (1/sqrt(2),0,1/sqrt(2))
				//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
				
				wert[0] += 2*1/sqrt(3)*3/2*h;
				wert[3] += 1/sqrt(3)*(-h)/2;
				wert[4] += 1/sqrt(3) *2*h;
				wert[5] += 1/sqrt(3)*(-h)/2;
				wert[6] += 1/sqrt(3) *2*h;
				
				//zeile[j] += 2*1/sqrt(3)*3/2*h;
				//zeile[j+1] += 1/sqrt(3)*(-h)/2;
				//zeile[j+2] += 1/sqrt(3) *2*h;
				//zeile[j+nx*ny] += 1/sqrt(3)*(-h)/2;
				//zeile[j+2*nx*ny] += 1/sqrt(3) *2*h;
			}
		
			else if(i % nx == 0) //rechte untere Kante ohne Eckpunkt
			{
				index[0]=j;
				index[1]=j+nx;
				index[2]=j-nx;
				index[3]=j-1;
				index[4]=j-2;
				index[5]=j+nx*ny;
				index[6]=j+2*nx*ny;
				
				//zentraler Differezenquotient in y-Richtung möglich
				wert[0] = -2;
				wert[1] = 1;
				wert[2] = 1;
				//modifizierter Differenzenquotient in x/z-Richtung
				wert[0] += 2*11/38;
				wert[3] = -28/38;
				wert[4] = 17/38;
				wert[5] = -28/38;
				wert[6] = 17/38;
				
				////zentraler Differezenquotient in y-Richtung möglich
				//zeile[j] = -2;
				//zeile[j+nx] = 1;
				//zeile[j-nx] = 1;
				////modifizierter Differenzenquotient in x/z-Richtung
				//zeile[j] += 2*11/38;
				//zeile[j-1] = -28/38;
				//zeile[j-2] = 17/38;
				//zeile[j+nx*ny] = -28/38;
				//zeile[j+2*nx*ny] = 17/38;
			
				//NeumannRB, Normalenvektor ist (-1/sqrt(2),0,1/sqrt(2))
				//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
				
				//wert[0] += 0*1/sqrt(3)*3/2*h;
				wert[3] += (-1)*1/sqrt(3)*(-h)/2;
				wert[4] += (-1)*1/sqrt(3) *2*h;
				wert[5] += 1/sqrt(3)*(-h)/2;
				wert[6] += 1/sqrt(3) *2*h;
				
				////zeile[j] += 0*1/sqrt(3)*3/2*h;
				//zeile[j-1] += (-1)*1/sqrt(3)*(-h)/2;
				//zeile[j-2] += (-1)*1/sqrt(3) *2*h;
				//zeile[j+nx*ny] += 1/sqrt(3)*(-h)/2;
				//zeile[j+2*nx*ny] += 1/sqrt(3) *2*h;
			}
		
			else // "innere" Punkte des Bodens
			{				
				index[0]=j;
				index[1]=j+1;
				index[2]=j-1;
				index[3]=j+nx;
				index[4]=j-nx;
				index[5]=j+nx*ny;
				index[6]=j+2*nx*ny;
				
				////zentraler Differenzenquotient in x/y-Richtung möglich
				
				wert[0] = -4;
				wert[1] = 1;
				wert[2] = 1;
				wert[3] = 1;
				wert[4] = 1;
				//modifizierter Differenzenquotient in z-Richtung
				wert[0] += 11/38;
				wert[5] = -28/38;
				wert[6] = 17/38;
				
				////zentraler Differenzenquotient in x/y-Richtung möglich
				//zeile[j] = -4;
				//zeile[j+1] = 1;
				//zeile[j-1] = 1;
				//zeile[j+nx] = 1;
				//zeile[j-nx] = 1;
				//modifizierter Differenzenquotient in z-Richtung
				//zeile[j] += 11/38;
				//zeile[j+nx*ny] = -28/38;
				//zeile[j+2*nx*ny] = 17/38;
			
				//NeumannRB, Normalenvektor ist (0,0,1))
				//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
				
				wert[0] += 3/2*h;
				wert[5] += (-h)/2;
				wert[6] += 2*h;
				
				//zeile[j] += 3/2*h;
				//zeile[j+nx*ny] += (-h)/2;
				//zeile[j+2*nx*ny] += 2*h;
		
			}
		}
	
		else if(i % nx*ny*nz > nx*ny*(nz-1)) //Deckel
		{
			if(i == 1) //vorderer oberer linker Eckpunkt
			{
				index[0]=j;
				index[1]=j+1;
				index[2]=j+2;
				index[3]=j+nx;
				index[4]=j+2*nx;
				index[5]=j-nx*ny;
				index[6]=j-2*nx*ny;
				
				//zentraler Differenzenquotient in keine Richtung möglich
				wert[0] = 3*11/38;
				wert[1] = -28/38;
				wert[2] = 17/38;
				wert[3] = -28/38;
				wert[4] = 17/38;
				wert[5] = -28/38;
				wert[6] = 17/38;
				
				////zentraler Differenzenquotient in keine Richtung möglich
				//zeile[j] = 3*11/38;
				//zeile[j+1] = -28/38;
				//zeile[j+2] = 17/38;
				//zeile[j+nx] = -28/38;
				//zeile[j+2*nx] = 17/38;
				//zeile[j-nx*ny] = -28/38;
				//zeile[j-2*nx*ny] = 17/38;
			
				//NeumannRB, Normalenvektor ist (1/sqrt(3),1/sqrt(3),-1/sqrt(3))
				//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
				
				wert[0] += 1*1/sqrt(3)*3/2*h;
				wert[1] += 1/sqrt(3)*(-h)/2;
				wert[2] += 1/sqrt(3) *2*h;
				wert[3] += 1/sqrt(3)*(-h)/2;
				wert[4] += 1/sqrt(3) *2*h;
				wert[5] += (-1)*1/sqrt(3)*(-h)/2;
				wert[6] += (-1)*1/sqrt(3) *2*h;
				
				//zeile[j] += 1*1/sqrt(3)*3/2*h;
				//zeile[j+1] += 1/sqrt(3)*(-h)/2;
				//zeile[j+2] += 1/sqrt(3) *2*h;
				//zeile[j+nx] += 1/sqrt(3)*(-h)/2;
				//zeile[j+2*nx] += 1/sqrt(3) *2*h;
				//zeile[j-nx*ny] += (-1)*1/sqrt(3)*(-h)/2;
				//zeile[j-2*nx*ny] += (-1)*1/sqrt(3) *2*h;
			}
		
			else if(i == nx) //vorderer oberer rechter Eckpunkt
			{
				index[0]=j;
				index[1]=j-1;
				index[2]=j-2;
				index[3]=j+nx;
				index[4]=j+2*nx;
				index[5]=j-nx*ny;
				index[6]=j-2*nx*ny;
				
				//zentraler Differenzenquotient in keine Richtung möglich
				wert[0] = 3*11/38;
				wert[1] = -28/38;
				wert[2] = 17/38;
				wert[3] = -28/38;
				wert[4] = 17/38;
				wert[5] = -28/38;
				wert[6] = 17/38;
				
				////zentraler Differenzenquotient in keine Richtung möglich
				//zeile[j] = 3*11/38;
				//zeile[j-1] = -28/38;
				//zeile[j-2] = 17/38;
				//zeile[j+nx] = -28/38;
				//zeile[j+2*nx] = 17/38;
				//zeile[j-nx*ny] = -28/38;
				//zeile[j-2*nx*ny] = 17/38;
			
				//NeumannRB, Normalenvektor ist (-1/sqrt(3),1/sqrt(3),-1/sqrt(3))
				//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
				
				wert[0] += (-1)*1/sqrt(3)*3/2*h;
				wert[1] += (-1)*1/sqrt(3)*(-h)/2;
				wert[2] += (-1)*1/sqrt(3) *2*h;
				wert[3] += 1/sqrt(3)*(-h)/2;
				wert[4] += 1/sqrt(3) *2*h;
				wert[5] += (-1)*1/sqrt(3)*(-h)/2;
				wert[6] += (-1)*1/sqrt(3) *2*h;
				
				//zeile[j] += (-1)*1/sqrt(3)*3/2*h;
				//zeile[j-1] += (-1)*1/sqrt(3)*(-h)/2;
				//zeile[j-2] += (-1)*1/sqrt(3) *2*h;
				//zeile[j+nx] += 1/sqrt(3)*(-h)/2;
				//zeile[j+2*nx] += 1/sqrt(3) *2*h;
				//zeile[j-nx*ny] += (-1)*1/sqrt(3)*(-h)/2;
				//zeile[j-2*nx*ny] += (-1)*1/sqrt(3) *2*h;
			}
		
			else if(i == nx*(ny-1)+1) //hinterer oberer linker Eckpunkt
			{
				index[0]=j;
				index[1]=j+1;
				index[2]=j+2;
				index[3]=j-nx;
				index[4]=j-2*nx;
				index[5]=j-nx*ny;
				index[6]=j-2*nx*ny;
				
				//zentraler Differenzenquotient in keine Richtung möglich
				wert[0] = 3*11/38;
				wert[1] = -28/38;
				wert[2] = 17/38;
				wert[3] = -28/38;
				wert[4] = 17/38;
				wert[5] = -28/38;
				wert[6] = 17/38;
				
				////zentraler Differenzenquotient in keine Richtung möglich
				//zeile[j] = 3*11/38;
				//zeile[j+1] = -28/38;
				//zeile[j+2] = 17/38;
				//zeile[j-nx] = -28/38;
				//zeile[j-2*nx] = 17/38;
				//zeile[j-nx*ny] = -28/38;
				//zeile[j-2*nx*ny] = 17/38;
			
				//NeumannRB, Normalenvektor ist (1/sqrt(3),-1/sqrt(3),-1/sqrt(3))
				//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
				
				wert[0] += (-1)*1/sqrt(3)*3/2*h;
				wert[1] += 1/sqrt(3)*(-h)/2;
				wert[2] += 1/sqrt(3) *2*h;
				wert[3] += (-1)*1/sqrt(3)*(-h)/2;
				wert[4] += (-1)*1/sqrt(3) *2*h;
				wert[5] += (-1)*1/sqrt(3)*(-h)/2;
				wert[6] += (-1)*1/sqrt(3) *2*h;
				
				//zeile[j] += (-1)*1/sqrt(3)*3/2*h;
				//zeile[j+1] += 1/sqrt(3)*(-h)/2;
				//zeile[j+2] += 1/sqrt(3) *2*h;
				//zeile[j-nx] += (-1)*1/sqrt(3)*(-h)/2;
				//zeile[j-2*nx] += (-1)*1/sqrt(3) *2*h;
				//zeile[j-nx*ny] += (-1)*1/sqrt(3)*(-h)/2;
				//zeile[j-2*nx*ny] += (-1)*1/sqrt(3) *2*h;
			}
		
			else if(i == nx*ny) //hinterer  oberer rechter Eckpunkt
			{
				index[0]=j;
				index[1]=j-1;
				index[2]=j-2;
				index[3]=j-nx;
				index[4]=j-2*nx;
				index[5]=j-nx*ny;
				index[6]=j-2*nx*ny;
				
				//zentraler Differenzenquotient in keine Richtung möglich
				wert[0] = 3*11/38;
				wert[1] = -28/38;
				wert[2] = 17/38;
				wert[3] = -28/38;
				wert[4] = 17/38;
				wert[5] = -28/38;
				wert[6] = 17/38;
				
				////zentraler Differenzenquotient in keine Richtung möglich
				//zeile[j] = 3*11/38;
				//zeile[j-1] = -28/38;
				//zeile[j-2] = 17/38;
				//zeile[j-nx] = -28/38;
				//zeile[j-2*nx] = 17/38;
				//zeile[j-nx*ny] = -28/38;
				//zeile[j-2*nx*ny] = 17/38;
			
				//NeumannRB, Normalenvektor ist (-1/sqrt(3),-1/sqrt(3),-1/sqrt(3))
				//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
				
				wert[0] += (-3)*1/sqrt(3)*3/2*h;
				wert[1] += (-1)*1/sqrt(3)*(-h)/2;
				wert[2] += (-1)*1/sqrt(3) *2*h;
				wert[3] += (-1)*1/sqrt(3)*(-h)/2;
				wert[4] += (-1)*1/sqrt(3) *2*h;
				wert[5] += (-1)*1/sqrt(3)*(-h)/2;
				wert[6] += (-1)*1/sqrt(3) *2*h;
				
				//zeile[j] += (-3)*1/sqrt(3)*3/2*h;
				//zeile[j-1] += (-1)*1/sqrt(3)*(-h)/2;
				//zeile[j-2] += (-1)*1/sqrt(3) *2*h;
				//zeile[j-nx] += (-1)*1/sqrt(3)*(-h)/2;
				//zeile[j-2*nx] += (-1)*1/sqrt(3) *2*h;
				//zeile[j-nx*ny] += (-1)*1/sqrt(3)*(-h)/2;
				//zeile[j-2*nx*ny] += (-1)*1/sqrt(3) *2*h;
			}
		
			else if(i % nx <= nx) //vordere obere Kante ohne Eckpunkt
			{
				index[0]=j;
				index[1]=j+1;
				index[2]=j-1;
				index[3]=j+nx;
				index[4]=j+2*nx;
				index[5]=j-nx*ny;
				index[6]=j-2*nx*ny;
				
				//zentraler Differezenquotient in x-Richtung möglich
				wert[0] = -2;
				wert[1] = 1;
				wert[2] = 1;
				//modifizierter Differenzenquotient in y/z-Richtung
				wert[0] += 2*11/38;
				wert[3] = -28/38;
				wert[4] = 17/38;
				wert[5] = -28/38;
				wert[6] = 17/38;
				
				////zentraler Differezenquotient in x-Richtung möglich
				//zeile[j] = -2;
				//zeile[j+1] = 1;
				//zeile[j-1] = 1;
				////modifizierter Differenzenquotient in y/z-Richtung
				//zeile[j] += 2*11/38;
				//zeile[j+nx] = -28/38;
				//zeile[j+2*nx] = 17/38;
				//zeile[j-nx*ny] = -28/38;
				//zeile[j-2*nx*ny] = 17/38;
			
				//NeumannRB, Normalenvektor ist (0,1/sqrt(2),-1/sqrt(2))
				//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
				
				//wert[0] += 0*1/sqrt(3)*3/2*h;
				wert[3] += 1/sqrt(3)*(-h)/2;
				wert[4] += 1/sqrt(3) *2*h;
				wert[5] += (-1)*1/sqrt(3)*(-h)/2;
				wert[6] += (-1)*1/sqrt(3) *2*h;
				
				////zeile[j] += 0*1/sqrt(3)*3/2*h;
				//zeile[j+nx] += 1/sqrt(3)*(-h)/2;
				//zeile[j+2*nx] += 1/sqrt(3) *2*h;
				//zeile[j-nx*ny] += (-1)*1/sqrt(3)*(-h)/2;
				//zeile[j-2*nx*ny] += (-1)*1/sqrt(3) *2*h;
			}
		
			else if(i % nx*ny > nx*(ny-1)) //hintere obere Kante ohne Eckpunkt
			{
				index[0]=j;
				index[1]=j+1;
				index[2]=j-1;
				index[3]=j-nx;
				index[4]=j-2*nx;
				index[5]=j-nx*ny;
				index[6]=j-2*nx*ny;
				
				//zentraler Differezenquotient in x-Richtung möglich
				wert[0] = -2;
				wert[1] = 1;
				wert[2] = 1;
				//modifizierter Differenzenquotient in y/z-Richtung
				wert[0] += 2*11/38;
				wert[3] = -28/38;
				wert[4] = 17/38;
				wert[5] = -28/38;
				wert[6] = 17/38;
				
				////zentraler Differezenquotient in x-Richtung möglich
				//zeile[j] = -2;
				//zeile[j+1] = 1;
				//zeile[j-1] = 1;
				////modifizierter Differenzenquotient in y/z-Richtung
				//zeile[j] += 2*11/38;
				//zeile[j-nx] = -28/38;
				//zeile[j-2*nx] = 17/38;
				//zeile[j-nx*ny] = -28/38;
				//zeile[j-2*nx*ny] = 17/38;
			
				//NeumannRB, Normalenvektor ist (0,-1/sqrt(2),-1/sqrt(2))
				//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
				
				wert[0] += (-2)*1/sqrt(3)*3/2*h;
				wert[3] += (-1)*1/sqrt(3)*(-h)/2;
				wert[4] += (-1)*1/sqrt(3) *2*h;
				wert[5] += (-1)*1/sqrt(3)*(-h)/2;
				wert[6] += (-1)*1/sqrt(3) *2*h;
				
				//zeile[j] += (-2)*1/sqrt(3)*3/2*h;
				//zeile[j-nx] += (-1)*1/sqrt(3)*(-h)/2;
				//zeile[j-2*nx] += (-1)*1/sqrt(3) *2*h;
				//zeile[j-nx*ny] += (-1)*1/sqrt(3)*(-h)/2;
				//zeile[j-2*nx*ny] += (-1)*1/sqrt(3) *2*h;
			}
		
			else if(i % nx == 1) //linke obere Kante ohne Eckpunkt
			{
				index[0]=j;
				index[1]=j+nx;
				index[2]=j-nx;
				index[3]=j+1;
				index[4]=j+2;
				index[5]=j-nx*ny;
				index[6]=j-2*nx*ny;
				
				//zentraler Differezenquotient in y-Richtung möglich
				wert[0] = -2;
				wert[1] = 1;
				wert[2] = 1;
				//modifizierter Differenzenquotient in x/z-Richtung
				wert[0] += 2*11/38;
				wert[3] = -28/38;
				wert[4] = 17/38;
				wert[5] = -28/38;
				wert[6] = 17/38;
				
				////zentraler Differezenquotient in y-Richtung möglich
				//zeile[j] = -2;
				//zeile[j+nx] = 1;
				//zeile[j-nx] = 1;
				////modifizierter Differenzenquotient in x/z-Richtung
				//zeile[j] += 2*11/38;
				//zeile[j+1] = -28/38;
				//zeile[j+2] = 17/38;
				//zeile[j-nx*ny] = -28/38;
				//zeile[j-2*nx*ny] = 17/38;
			
				//NeumannRB, Normalenvektor ist (1/sqrt(2),0,-1/sqrt(2))
				//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
				
				//wert[0] += 0*1/sqrt(3)*3/2*h;
				wert[3] += 1/sqrt(3)*(-h)/2;
				wert[4] += 1/sqrt(3) *2*h;
				wert[5] += (-1)*1/sqrt(3)*(-h)/2;
				wert[6] += (-1)*1/sqrt(3) *2*h;
				
				////zeile[j] += 0*1/sqrt(3)*3/2*h;
				//zeile[j+1] += 1/sqrt(3)*(-h)/2;
				//zeile[j+2] += 1/sqrt(3) *2*h;
				//zeile[j-nx*ny] += (-1)*1/sqrt(3)*(-h)/2;
				//zeile[j-2*nx*ny] += (-1)*1/sqrt(3) *2*h;
			}
		
			else if(i % nx == 0) //rechte obere Kante ohne Eckpunkt
			{
				index[0]=j;
				index[1]=j+nx;
				index[2]=j-nx;
				index[3]=j-1;
				index[4]=j-2;
				index[5]=j-nx*ny;
				index[6]=j-2*nx*ny;
				
				//zentraler Differezenquotient in y-Richtung möglich
				wert[0] = -2;
				wert[1] = 1;
				wert[2] = 1;
				//modifizierter Differenzenquotient in x/z-Richtung
				wert[0] += 2*11/38;
				wert[3] = -28/38;
				wert[4] = 17/38;
				wert[5] = -28/38;
				wert[6] = 17/38;
				
				////zentraler Differezenquotient in y-Richtung möglich
				//zeile[j] = -2;
				//zeile[j+nx] = 1;
				//zeile[j-nx] = 1;
				////modifizierter Differenzenquotient in x/z-Richtung
				//zeile[j] += 2*11/38;
				//zeile[j-1] = -28/38;
				//zeile[j-2] = 17/38;
				//zeile[j-nx*ny] = -28/38;
				//zeile[j-2*nx*ny] = 17/38;
			
				//NeumannRB, Normalenvektor ist (-1/sqrt(2),0,-1/sqrt(2))
				//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
				
				wert[0] += (-2)*1/sqrt(3)*3/2*h;
				wert[3] += (-1)*1/sqrt(3)*(-h)/2;
				wert[4] += (-1)*1/sqrt(3) *2*h;
				wert[5] += (-1)*1/sqrt(3)*(-h)/2;
				wert[6] += (-1)*1/sqrt(3) *2*h;
				
				//zeile[j] += (-2)*1/sqrt(3)*3/2*h;
				//zeile[j-1] += (-1)*1/sqrt(3)*(-h)/2;
				//zeile[j-2] += (-1)*1/sqrt(3) *2*h;
				//zeile[j-nx*ny] += (-1)*1/sqrt(3)*(-h)/2;
				//zeile[j-2*nx*ny] += (-1)*1/sqrt(3) *2*h;
			}
		
			else // "innere" Punkte des Deckels
			{
				index[0]=j;
				index[1]=j+1;
				index[2]=j-1;
				index[3]=j+nx;
				index[4]=j-nx;
				index[5]=j-nx*ny;
				index[6]=j-2*nx*ny;
				
				//zentraler Differenzenquotient in x/y-Richtung möglich
				wert[0] = -4;
				wert[1] = 1;
				wert[2] = 1;
				wert[3] = 1;
				wert[4] = 1;
				//modifizierter Differenzenquotient in z-Richtung
				wert[0] += 11/38;
				wert[5] = -28/38;
				wert[6] = 17/38;
				
				////zentraler Differenzenquotient in x/y-Richtung möglich
				//zeile[j] = -4;
				//zeile[j+1] = 1;
				//zeile[j-1] = 1;
				//zeile[j+nx] = 1;
				//zeile[j-nx] = 1;
				////modifizierter Differenzenquotient in z-Richtung
				//zeile[j] += 11/38;
				//zeile[j-nx*ny] = -28/38;
				//zeile[j-2*nx*ny] = 17/38;
			
				//NeumannRB, Normalenvektor ist (0,0,-1))
				//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
				
				wert[0] += (-1)*3/2*h;
				wert[5] += (-1)*(-h)/2;
				wert[6] += (-1)*2*h;
				
				//zeile[j] += (-1)*3/2*h;
				//zeile[j-nx*ny] += (-1)*(-h)/2;
				//zeile[j-2*nx*ny] += (-1)*2*h;
		
			}
		}	
    
		else if(i % nx*ny <= nx) //vordere Seite, aber nicht Boden oder Deckel
		{
			if(i % nx==1) //linke Kante
			{
				index[0]=j;
				index[1]=j+nx*ny;
				index[2]=j-nx*ny;
				index[3]=j+1;
				index[4]=j+2;
				index[5]=j+nx;
				index[6]=j+2*nx;
				
				//zentraler Differenzenquotient nur in z-Richtung möglich
				wert[0] = -2;
				wert[1] = 1;
				wert[2] = 1;
				//modifizierter Differenzenquotient in x/y-Richtung
				wert[0] += 2*11/38;
				wert[3] = -28/38;
				wert[4] = 17/38;
				wert[5] = -28/38;
				wert[6] = 17/38;
				
				////zentraler Differenzenquotient nur in z-Richtung möglich
				//zeile[j] = -2;
				//zeile[j+nx*ny] = 1;
				//zeile[j-nx*ny] = 1;
				////modifizierter Differenzenquotient in x/y-Richtung
				//zeile[j] += 2*11/38;
				//zeile[j+1] = -28/38;
				//zeile[j+2] = 17/38;
				//zeile[j+nx] = -28/38;
				//zeile[j+2*nx] = 17/38;
			
				//NeumannRB, Normalenvektor ist (1/sqrt(2),1/sqrt(2),0)
				//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
				
				wert[0] += 2*1/sqrt(3)*3/2*h;
				wert[3] += 1/sqrt(3)*(-h)/2;
				wert[4] += 1/sqrt(3) *2*h;
				wert[5] += 1/sqrt(3)*(-h)/2;
				wert[6] += 1/sqrt(3) *2*h;
				
				//zeile[j] += 2*1/sqrt(3)*3/2*h;
				//zeile[j+1] += 1/sqrt(3)*(-h)/2;
				//zeile[j+2] += 1/sqrt(3) *2*h;
				//zeile[j+nx] += 1/sqrt(3)*(-h)/2;
				//zeile[j+2*nx] += 1/sqrt(3) *2*h;
			
			}

			if(i % nx==0) //rechte Kante
			{
				index[0]=j;
				index[1]=j+nx*ny;
				index[2]=j-nx*ny;
				index[3]=j-1;
				index[4]=j-2;
				index[5]=j+nx;
				index[6]=j+2*nx;
				
				//zentraler Differenzenquotient nur in z-Richtung möglich
				wert[0] = -2;
				wert[1] = 1;
				wert[2] = 1;
				//modifizierter Differenzenquotient in x/y-Richtung
				wert[0] += 2*11/38;
				wert[3] = -28/38;
				wert[4] = 17/38;
				wert[5] = -28/38;
				wert[6] = 17/38;
				
				////zentraler Differenzenquotient nur in z-Richtung möglich
				//zeile[j] = -2;
				//zeile[j+nx*ny] = 1;
				//zeile[j-nx*ny] = 1;
				////modifizierter Differenzenquotient in x/y-Richtung
				//zeile[j] += 2*11/38;
				//zeile[j-1] = -28/38;
				//zeile[j-2] = 17/38;
				//zeile[j+nx] = -28/38;
				//zeile[j+2*nx] = 17/38;

				//NeumannRB, Normalenvektor ist (-1/sqrt(2),1/sqrt(2),0)
				//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
				
				//wert[0] += 0*1/sqrt(3)*3/2*h;
				wert[3] += (-1)*1/sqrt(3)*(-h)/2;
				wert[4] += (-1)*1/sqrt(3) *2*h;
				wert[5] += 1/sqrt(3)*(-h)/2;
				wert[6] += 1/sqrt(3) *2*h;
				
				////zeile[j] += 0*1/sqrt(3)*3/2*h;
				//zeile[j-1] += (-1)*1/sqrt(3)*(-h)/2;
				//zeile[j-2] += (-1)*1/sqrt(3) *2*h;
				//zeile[j+nx] += 1/sqrt(3)*(-h)/2;
				//zeile[j+2*nx] += 1/sqrt(3) *2*h;			
			}
		
			else //vordere "innere" Seite
			{
				index[0]=j;
				index[1]=j+nx*ny;
				index[2]=j-nx*ny;
				index[3]=j+1;
				index[4]=j-1;
				index[5]=j+nx;
				index[6]=j+2*nx;
				
				//zentraler Differenzenquotient in x/z-Richtung möglich
				wert[0] = -4;
				wert[1] = 1;
				wert[2] = 1;
				wert[3]=1;
				wert[4]=1;
				//modifizierter Differenzenquotient in y-Richtung
				wert[0] += 11/38;
				wert[5] = -28/38;
				wert[6] = 17/38;
								
				////zentraler Differenzenquotient in x/z-Richtung möglich
				//zeile[j] = -4;
				//zeile[j+nx*ny] = 1;
				//zeile[j-nx*ny] = 1;
				//zeile[j+1]=1;
				//zeile[j-1]=1;
				////modifizierter Differenzenquotient in y-Richtung
				//zeile[j] += 11/38;
				//zeile[j+nx] = -28/38;
				//zeile[j+2*nx] = 17/38;
			
				//NeumannRB, Normalenvektor ist (0,1,0)
				//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
				
				wert[0] += 3/2*h;
				wert[5] += (-h)/2;
				wert[6] += 2*h;
				
				//zeile[j] += 3/2*h;
				//zeile[j+nx] += (-h)/2;
				//zeile[j+2*nx] += 2*h;
			}
		}
	
		else if(i % nx*ny > nx*(ny-1)) //hintere Seite, aber nicht Boden oder Deckel
		{
			if(i % nx==1) //linke Kante
			{
				index[0]=j;
				index[1]=j+nx*ny;
				index[2]=j-nx*ny;
				index[3]=j+1;
				index[4]=j+2;
				index[5]=j-nx;
				index[6]=j-2*nx;
				
				//zentraler Differenzenquotient nur in z-Richtung möglich
				wert[0] = -2;
				wert[1] = 1;
				wert[2] = 1;
				//modifizierter Differenzenquotient in x/y-Richtung
				wert[0] += 2*11/38;
				wert[3] = -28/38;
				wert[4] = 17/38;
				wert[5] = -28/38;
				wert[6] = 17/38;
				
				////zentraler Differenzenquotient nur in z-Richtung möglich
				//zeile[j] = -2;
				//zeile[j+nx*ny] = 1;
				//zeile[j-nx*ny] = 1;
				////modifizierter Differenzenquotient in x/y-Richtung
				//zeile[j] += 2*11/38;
				//zeile[j+1] = -28/38;
				//zeile[j+2] = 17/38;
				//zeile[j-nx] = -28/38;
				//zeile[j-2*nx] = 17/38;
			
				//NeumannRB, Normalenvektor ist (1/sqrt(2),-1/sqrt(2),0)
				//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
				
				//wert[0] += 0*1/sqrt(3)*3/2*h;
				wert[3] += 1/sqrt(3)*(-h)/2;
				wert[4] += 1/sqrt(3) *2*h;
				wert[5] += (-1)*1/sqrt(3)*(-h)/2;
				wert[6] += (-1)*1/sqrt(3) *2*h;
				
				////zeile[j] += 0*1/sqrt(3)*3/2*h;
				//zeile[j+1] += 1/sqrt(3)*(-h)/2;
				//zeile[j+2] += 1/sqrt(3) *2*h;
				//zeile[j-nx] += (-1)*1/sqrt(3)*(-h)/2;
				//zeile[j-2*nx] += (-1)*1/sqrt(3) *2*h;
			
			}

			if(i % nx==0) //rechte Kante
			{
				index[0]=j;
				index[1]=j+nx*ny;
				index[2]=j-nx*ny;
				index[3]=j-1;
				index[4]=j-2;
				index[5]=j-nx;
				index[6]=j-2*nx;
				
				//zentraler Differenzenquotient nur in z-Richtung möglich
				wert[0] = -2;
				wert[1] = 1;
				wert[2] = 1;
				//modifizierter Differenzenquotient in x/y-Richtung
				wert[0] += 2*11/38;
				wert[3] = -28/38;
				wert[4] = 17/38;
				wert[5] = -28/38;
				wert[6] = 17/38;
				
				////zentraler Differenzenquotient nur in z-Richtung möglich
				//zeile[j] = -2;
				//zeile[j+nx*ny] = 1;
				//zeile[j-nx*ny] = 1;
				////modifizierter Differenzenquotient in x/y-Richtung
				//zeile[j] += 2*11/38;
				//zeile[j-1] = -28/38;
				//zeile[j-2] = 17/38;
				//zeile[j-nx] = -28/38;
				//zeile[j-2*nx] = 17/38;
			
				//NeumannRB, Normalenvektor ist (-1/sqrt(2),-1/sqrt(2),0)
				//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
				wert[0] += (-2)*1/sqrt(3)*3/2*h;
				wert[3] += (-1)*1/sqrt(3)*(-h)/2;
				wert[4] += (-1)*1/sqrt(3) *2*h;
				wert[5] += (-1)*1/sqrt(3)*(-h)/2;
				wert[6] += (-1)*1/sqrt(3) *2*h;
				
				//zeile[j] += (-2)*1/sqrt(3)*3/2*h;
				//zeile[j-1] += (-1)*1/sqrt(3)*(-h)/2;
				//zeile[j-2] += (-1)*1/sqrt(3) *2*h;
				//zeile[j-nx] += (-1)*1/sqrt(3)*(-h)/2;
				//zeile[j-2*nx] += (-1)*1/sqrt(3) *2*h;
			}
		
			else // hintere "innere" Seite
			{
				index[0]=j;
				index[1]=j+nx*ny;
				index[2]=j-nx*ny;
				index[3]=j+1;
				index[4]=j-1;
				index[5]=j-nx;
				index[6]=j-2*nx;
				
				//zentraler Differenzenquotient in x/z-Richtung möglich
				wert[0] = -4;
				wert[1] = 1;
				wert[2] = 1;
				wert[3] =1;
				wert[4] =1;
				//modifizierter Differenzenquotient in y-Richtung
				wert[0] += 11/38;
				wert[5] = -28/38;
				wert[6] = 17/38;
				
				////zentraler Differenzenquotient in x/z-Richtung möglich
				//zeile[j] = -4;
				//zeile[j+nx*ny] = 1;
				//zeile[j-nx*ny] = 1;
				//zeile[j+1] =1;
				//zeile[j-1] =1;
				////modifizierter Differenzenquotient in y-Richtung
				//zeile[j] += 11/38;
				//zeile[j-nx] = -28/38;
				//zeile[j-2*nx] = 17/38;
				
				//NeumannRB, Normalenvektor ist (0,-1,0)
				//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
				wert[0] += (-1)*3/2*h;
				wert[5] += (-1)*(-h)/2;
				wert[6] += (-1)*2*h;
				
				//zeile[j] += (-1)*3/2*h;
				//zeile[j-nx] += (-1)*(-h)/2;
				//zeile[j-2*nx] += (-1)*2*h;
			}	
		}	
	
		else if(i % nx ==1) //linke Seite, aber nicht vordere/hintere Seite oder Boden/Deckel
		{
			index[0]=j;
			index[1]=j+nx;
			index[2]=j-nx;
			index[3]=j+nx*ny;
			index[4]=j-nx*ny;
			index[5]=j+1;
			index[6]=j+2;
			
			//zentraler Differenzenquotient ist nur in x-Richtung nicht möglich, deshalb zuerst normal in y/z-Richtung
			wert[0] = -4;
			wert[1] = 1;
			wert[2] = 1;
			wert[3] = 1;
			wert[4] = 1;
			//modifizierter Differenzenquotient in x-Richtung
			wert[0] += 11/38;
			wert[5] = -28/38;
			wert[6] = 17/38;
			
			////zentraler Differenzenquotient ist nur in x-Richtung nicht möglich, deshalb zuerst normal in y/z-Richtung
			//zeile[j] = -4;
			//zeile[j+nx] = 1;
			//zeile[j-nx] = 1;
			//zeile[j+nx*ny] = 1;
			//zeile[j-nx*ny] = 1;
			////modifizierter Differenzenquotient in x-Richtung
			//zeile[j] += 11/38;
			//zeile[j+1] = -28/38;
			//zeile[j+2] = 17/38;
		
			//NeumannRB, Normalenvektor ist (1,0,0)
			//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
			
			wert[0] += 3/2*h;
			wert[5] += (-h)/2;
			wert[6] += 2*h;
			
			//zeile[j] += 3/2*h;
			//zeile[j+1] += (-h)/2;
			//zeile[j+2] += 2*h;
		}
		
		else if(i % nx ==0) //rechte Seite, aber nicht vordere/hintere Seite oder Boden/Deckel
		{
			index[0]=j;
			index[1]=j+nx;
			index[2]=j-nx;
			index[3]=j+nx*ny;
			index[4]=j-nx*ny;
			index[5]=j-1;
			index[6]=j-2;
			
			//zentraler Differenzenquotient ist nur in x-Richtung nicht möglich, deshalb zuerst normal in y/z-Richtung
			wert[0] = -4;
			wert[1] = 1;
			wert[2] = 1;
			wert[3] = 1;
			wert[4] = 1;
			//modifizierter Differenzenquotient in x-Richtung
			wert[0] += 11/38;
			wert[5] = -28/38;
			wert[6] = 17/38;
						
			////zentraler Differenzenquotient ist nur in x-Richtung nicht möglich, deshalb zuerst normal in y/z-Richtung
			//zeile[j] = -4;
			//zeile[j+nx] = 1;
			//zeile[j-nx] = 1;
			//zeile[j+nx*ny] = 1;
			//zeile[j-nx*ny] = 1;
			////modifizierter Differenzenquotient in x-Richtung
			//zeile[j] += 11/38;
			//zeile[j-1] = -28/38;
			//zeile[j-2] = 17/38;
		
			//NeumannRB, Normalenvektor ist (-1,0,0)
			//RB wird auf die normale Zeile addiert, um die quadratische Struktur beizubehalten
			
			wert[0] += (-1)*3/2*h;
			wert[5] += (-1)*(-h)/2;
			wert[6] += (-1)*2*h;
			
			//zeile[j] += (-1)*3/2*h;
			//zeile[j-1] += (-1)*(-h)/2;
			//zeile[j-2] += (-1)*2*h;
		}
	
		else //innere Punkte
		{    
			index[0]=j;
			index[1]=j+1;
			index[2]=j-1;
			index[3]=j+nx;
			index[4]=j-nx;
			index[5]=j+nx*ny;
			index[6]=j-nx*ny;
			
			//zentraler Differenzenquotient in alle Richtung möglich
			wert[0] = -6;
			wert[1] = 1;
			wert[2] = 1;
			wert[3] = 1;
			wert[4] = 1;
			wert[5] = 1;
			wert[6] = 1;
			
			////zentraler Differenzenquotient in alle Richtung möglich
			//zeile[j] = -6;
			//zeile[j+1] = 1;
			//zeile[j-1] = 1;
			//zeile[j+nx] = 1;
			//zeile[j-nx] = 1;
			//zeile[j+nx*ny] = 1;
			//zeile[j-nx*ny] = 1;
			
			//keine RB
		}	
		
		A.sequential_fill(index[0],wert[0]);
		A.sequential_fill(index[1],wert[1]);
		A.sequential_fill(index[2],wert[2]);
		A.sequential_fill(index[3],wert[3]);
		A.sequential_fill(index[4],wert[4]);
		A.sequential_fill(index[5],wert[5]);
		A.sequential_fill(index[6],wert[6]);
		
		A.end_of_row();
	}
}

}
#endif // __ASSEMBLE_HPP_
