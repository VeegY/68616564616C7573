// #################################################################################################
//			Studienprojekt Modellbildung & Simulation - 2015/16
// #################################################################################################
// 					Header: gpumatoperations.hpp
// ------------------------------------Doxygen-Dokumentation----------------------------------------
///  \file gpumatoperations.hpp
///  \brief
///  Stellt Matrix Vektor Operationen auf der GPU zur Verfügung.
///  Funktionen müssen mit nvcc Compiler kompiliert werden.
//#################################################################################################
#	pragma once
#	include <iostream>
#	include <cstdlib>
#	include <ctime>
#	include <cuda_runtime.h>
#   include "DIA.hpp"


//kernel für spmv(): Ein Thread pro Reihe
template <typename datatype>
__global__ void spmvkernel(datatype* res, datatype* matval, int* matoff, datatype* rhs, int ndiags, int n)
{
	int idx = blockIdx.x*blockDim.x + threadIdx.x;
	int colind;
	if (idx < n){
		res[idx] = 0;
		for (int i = 0; i < ndiags; i++)
		{
			colind = idx + matoff[i];
			if ( colind >= 0 && colind < n){
				res[idx] += matval[idx + i*n] * rhs[colind];
			}
		}
	}
}

//kernel für gpudefect(): Ein Thread pro Reihe
template <typename datatype>
__global__ void defectkernel(datatype* res, datatype* matval, int* matoff,datatype* vec, datatype* rhs, int ndiags, int n)
{
	int idx = blockIdx.x*blockDim.x + threadIdx.x;
	int colind;
	if (idx < n){
		res[idx] = vec[idx];
		for (int i = 0; i < ndiags; i++)
		{
			colind = idx + matoff[i];
			if (colind >= 0 && colind < n){
				res[idx] -= matval[idx + i*n] * rhs[colind];
			}
		}
	}
}

// Sparse Matrix Vektor Multiplikation
// Bsp: spmv(r, mat, rhs);  brechnet r=mat*rhs
// Dimension der Matrix ist dabei auf 65535*512 begrenzt (max Blocks per Dimension * max Threads per Blocks)
template <typename datatype>
void spmv(Vector<datatype>& res, DIA<datatype>& mat, Vector<datatype>& rhs){
	cudaDeviceProp devProp;
	cudaGetDeviceProperties(&devProp, 0);
	int maxthread = devProp.maxThreadsPerBlock;
	int dim = mat.dim();
	int *d_matoff;
	datatype *d_res;
	datatype *d_matval;
	datatype *d_rhs;
	datatype *h_solution = new datatype[dim];
	if (res._dim != mat._dim || mat._dim != rhs._dim)
	{
		throw invalid_argument(" -Achtung! Dimensionsfehler!- ");
	}
	else{
		/*========Speicherallokation=========*/
		if (cudaSuccess != cudaMalloc(&d_matval, sizeof(datatype)*dim*mat._numDiags))
		{
			cout << "allocate error" << endl;

		}

		if (cudaSuccess != cudaMalloc(&d_rhs, sizeof(datatype)*dim))
		{
			cout << "allocate error" << endl;

		}

		if (cudaSuccess != cudaMalloc(&d_res, sizeof(datatype)*dim))
		{
			cout << "allocate error" << endl;
		}

		if (cudaSuccess != cudaMalloc(&d_matoff, sizeof(int)*mat._numDiags))
		{
			cout << "allocate error" << endl;
		}
		/*========Memcpy=========*/
		if (cudaSuccess != cudaMemcpy(d_matval, mat._data->_data, sizeof(datatype)*dim*mat._numDiags, cudaMemcpyHostToDevice))
		{
			cout << "failed to copy" << endl;
		}

		if (cudaSuccess != cudaMemcpy(d_matoff, mat._offset->_data, sizeof(int)*mat._numDiags, cudaMemcpyHostToDevice))
		{
			cout << "failed to copy" << endl;
		}

		if (cudaSuccess != cudaMemcpy(d_res, res._data, sizeof(datatype)*dim, cudaMemcpyHostToDevice))
		{
			cout << "failed to copy" << endl;
		}

		if (cudaSuccess != cudaMemcpy(d_rhs, rhs._data, sizeof(datatype)*dim, cudaMemcpyHostToDevice))
		{
			cout << "failed to copy" << endl;
		}
		//Kernel-Launch
		spmvkernel<<<(dim/maxthread)+1, maxthread>>>(d_res, d_matval, d_matoff, d_rhs, mat._numDiags, dim);


		if (cudaSuccess != cudaGetLastError())
		{
			cout << "kernel launch failed" << endl;
		}

		//Ergebnis kopieren
		if (cudaSuccess != cudaMemcpy(res._data, d_res, sizeof(datatype)*dim, cudaMemcpyDeviceToHost))
		{
			cout << "failed to copy" << endl;
		}
		//Free
		cudaFree(d_matval);
		cudaFree(d_matoff);
		cudaFree(d_res);
		cudaFree(d_rhs);

	}
}


// Defektberechnung
// Bsp: gpudefect(res, mat, rhs, vec); setzt res = vec - mat*rhs
template <typename datatype>
void gpudefect(Vector<datatype>& res, DIA<datatype>& mat, Vector<datatype>& vec, Vector<datatype>& rhs){
	int dim = mat.dim();
	int *d_matoff;
	cudaDeviceProp devProp;
	cudaGetDeviceProperties(&devProp, 0);
	int maxthread = devProp.maxThreadsPerBlock;
	datatype *d_res;
	datatype *d_matval;
	datatype *d_rhs;
	datatype *d_vec;
	datatype *h_solution = new datatype[dim];
	if (res._dim != mat._dim || mat._dim != rhs._dim)
	{
		throw invalid_argument(" -Achtung! Dimensionsfehler!- ");
	}
	else{
		/*========Speicherallokation=========*/
		if (cudaSuccess != cudaMalloc(&d_matval, sizeof(datatype)*dim*mat._numDiags))
		{
			cout << "allocate error" << endl;
		}

		if (cudaSuccess != cudaMalloc(&d_rhs, sizeof(datatype)*dim))
		{
			cout << "allocate error" << endl;
		}

		if (cudaSuccess != cudaMalloc(&d_res, sizeof(datatype)*dim))
		{
			cout << "allocate error" << endl;
		}

		if (cudaSuccess != cudaMalloc(&d_vec, sizeof(datatype)*dim))
		{
			cout << "allocate error" << endl;
		}

		if (cudaSuccess != cudaMalloc(&d_matoff, sizeof(int)*mat._numDiags))
		{
			cout << "allocate error" << endl;
		}
		/*========Memcpy=========*/
		if (cudaSuccess != cudaMemcpy(d_matval, mat._data->_data, sizeof(datatype)*dim*mat._numDiags, cudaMemcpyHostToDevice))
		{
			cout << "failed to copy" << endl;
		}

		if (cudaSuccess != cudaMemcpy(d_matoff, mat._offset->_data, sizeof(int)*mat._numDiags, cudaMemcpyHostToDevice))
		{
			cout << "failed to copy" << endl;
		}

		if (cudaSuccess != cudaMemcpy(d_res, res._data, sizeof(datatype)*dim, cudaMemcpyHostToDevice))
		{
			cout << "failed to copy" << endl;
		}

		if (cudaSuccess != cudaMemcpy(d_rhs, rhs._data, sizeof(datatype)*dim, cudaMemcpyHostToDevice))
		{
			cout << "failed to copy" << endl;
		}

		if (cudaSuccess != cudaMemcpy(d_vec, vec._data, sizeof(datatype)*dim, cudaMemcpyHostToDevice))
		{
			cout << "failed to copy" << endl;
		}
		//Kernel-Launch
		defectkernel << <(dim / maxthread) + 1, maxthread >> >(d_res, d_matval, d_matoff, d_vec, d_rhs, mat._numDiags, dim);


		if (cudaSuccess != cudaGetLastError())
		{
			cout << "kernel launch failed" << endl;
		}

		//Ergebnis kopieren
		if (cudaSuccess != cudaMemcpy(res._data, d_res, sizeof(datatype)*dim, cudaMemcpyDeviceToHost))
		{
			cout << "failed to copy" << endl;
		}
		//Free
		cudaFree(d_matval);
		cudaFree(d_matoff);
		cudaFree(d_res);
		cudaFree(d_rhs);
		cudaFree(d_vec);
	}
}
