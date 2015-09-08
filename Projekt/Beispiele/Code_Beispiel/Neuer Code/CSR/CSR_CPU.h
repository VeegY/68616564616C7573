// ################################################################################################
//					 FD - Code Project
// ################################################################################################
//					  Header: CSR_CPU
// ------------------------------------Doxygen-Dokumentation---------------------------------------
///  \file CSR_CPU.h
///  \brief
///  CPU Implementierung der CSR-Matrix-Formats unter Verwendung der 'vector'-Implementierung.
// ------------------------------------------------------------------------------------------------
// Verwendete Header:
#    pragma once
#    include <iostream>

#    include "Miscellaneous.h"
#    include "vector_CPU.h"
// -------------------------------------------------------------------------------------------------
// Verwendet in :
//   CSR_CPU_TEST.h
// -------------------------------------------------------------------------------------------------
// Copyright : 2015 malte-christian.loddoch@tu-dortmund.de
// #################################################################################################

// ========================================DOKUMENTATION============================================
///\class 	csr                         
///\brief 	Implementiert das CSR-Matrixformat.
// -------------------------------------------------------------------------------------------------
///@tparam	dtype 	Template Type fuer die 'csr'-Klasse
// =================================================================================================
   template<typename dtype> ///
   class csr
   {
public:

#pragma region PUBLIC_VARIABLES
	int _dim;	///<	_dim 		Dimension der Matrix.
	int _nnvals;	///<	_nnvals 	Anzahl der Nichtnulleintraege.
	dtype* _data;	///<	_data		Array der Matrix-Nichtnulleintraege.
	int* _col_idx;	///<	_colidx		Array der Spaltenpositionen mit Laenge '_nnvals'.
	int* _row_idx;	///<	_rowidx		Array der Zeilenuebergang-Nichtnulleintraege mit Laenge '_dim + 1'.
#pragma endregion

#pragma region CONSTRUCTORS
// ========================================DOKUMENTATION============================================
//					Standard-Konstruktor
///\brief	Definition des Standard-Konstruktors.
// =================================================================================================
	csr() : _data(NULL), _col_idx(NULL), _row_idx(NULL) {}
// =================================================================================================

// ========================================DOKUMENTATION============================================
//					Direkt-Konstruktor
///\brief	Definition eines Konstruktors fuer direkte Zuweisungen
// -------------------------------------------------------------------------------------------------
///@param	dim 		Dimension der Matrix.
///@param	nnvals 		Anzahl der Nichtnulleintraege.
///@param	data		Array der Matrix-Nichtnulleintraege.
///@param	col_idx		Array der Spaltenpositionen mit Laenge '_nnvals'.
///@param	row_idx		Array der Zeilenuebergang-Nichtnulleintraege mit Laenge '_dim + 1'.
// -------------------------------------------------------------------------------------------------
	csr(int dim, int nnvals, dtype* data, int* col_idx, int* row_idx) : _dim(dim), _nnvals(nnvals),
					_data(data), _col_idx(col_idx), _row_idx(row_idx) {}
// =================================================================================================

// ========================================DOKUMENTATION============================================
//					   0-Konstruktor
///\brief	Definition eines Konstruktors fuer 0-Matrizen der Groesse 'dim'.
// -------------------------------------------------------------------------------------------------
///@param	dim 		Dimension der Matrix.
// -------------------------------------------------------------------------------------------------
	csr(int dim) : _dim(dim), _nnvals(0), _data(dtype(0.0)), _row_idx(0)
	{
		_col_idx = new int[dim];
		#pragma omp parallel for default(shared)
		for (int i = 0; i < dim + 1; i++) { _col_idx[i] = 0; }
	}
// =================================================================================================

// ========================================DOKUMENTATION============================================
//					 Copy-Konstruktor
///\brief	Definition des Copy-Konstruktors
// -------------------------------------------------------------------------------------------------
///@param	init_csr	Referenz der zu kopierenden Matriz.
// -------------------------------------------------------------------------------------------------
	// const induziert einen Conversion-Fehler in ref_fex?! (aber nicht im =-op)
	csr(csr& init_csr) : _dim(init_csr._dim), _nnvals(init_csr._nnvals)
	{   
		// get proper pointers to the data, otherwise, compilers cannot vectorise the loop
		const dtype* __MYRESTRICT data = init_csr.ref_fed();
		const int* __MYRESTRICT col = init_csr.ref_fec();
		const int* __MYRESTRICT row = init_csr.ref_fer();

		// allocate data
		_data = new dtype[_nnvals];
		_col_idx = new int[_nnvals];
		_row_idx = new int[_dim + 1];

		// match data

		#pragma omp for schedule(static)
		for (int i = 0; i < _nnvals; i++)
		{
			_data[i] = data[i];
			_col_idx[i] = col[i];
		}
		#pragma omp for schedule(static)
		for (int i = 0; i < _dim + 1; i++) _row_idx[i] = row[i];
	}
// =================================================================================================

// ========================================DOKUMENTATION============================================
//					    Destruktor
///\brief	Definition des Destruktors
// -------------------------------------------------------------------------------------------------
	~csr()
	{ 
		// this invokes an assertion debug fail 
		//delete[] _data; 
		//delete[] _col_idx;
		//delete[] _row_idx;
	}
// =================================================================================================
#pragma endregion

#pragma region OPERATIONS
// ========================================DOKUMENTATION============================================
//					   '()'-Operator
///\brief	Ueberladung des '()'-Operator
// -------------------------------------------------------------------------------------------------
///@param	column		Array der Spaltenpositionen mit Laenge '_nnvals'.
///@param	row		Spaltenindex der Position des angeforderten Elements.
// -------------------------------------------------------------------------------------------------
	dtype operator() (const int column, const int row)
	{
		dtype result(dtype(0));

		//get indexes of columns for given row 
		int start_id(_row_idx[row]);
		int end_id(_row_idx[row + 1]);

		//search all entries of '_col_idx' for 'column' in given row
		while (start_id < end_id)
		{
			if (_col_idx[start_id] == column) return _data[start_id];
			start_id++;
		}

		//return 0
		return result;
	}
// =================================================================================================

// ========================================DOKUMENTATION============================================
//				       Assignment-Operators
///\brief	Definition des Assignment-Operators
// -------------------------------------------------------------------------------------------------
///@param	mat 		Matrix die zugewiesen werden soll.
// -------------------------------------------------------------------------------------------------
	csr<dtype>& operator= (csr<dtype>& mat)
	{
		// check if Dimensions fit
		if (_dim == mat._dim)
		{
			// get proper pointers to the data, otherwise, compilers cannot vectorise the loop
			const dtype* __MYRESTRICT data = mat.ref_fed();
			const int* __MYRESTRICT col = mat.ref_fec();
			const int* __MYRESTRICT row = mat.ref_fer();

			// allocate memory
			_nnvals = mat._nnvals;
			_data = new dtype[_nnvals];
			_col_idx = new dtype[_nnvals];
			_row_idx = new dtype[_dim+1];

			// match data
			#pragma omp for schedule(static)
			for (int i = 0; i < _nnvals; i++)
			{
				_data[i] = data[i];
				_col_idx[i] = col[i];
			}
			#pragma omp for schedule(static)
			for (int i = 0; i < _dim + 1; i++) _row_idx[i] = row[i];
			return *this;
		}
		else
		{
			std::cout << "Dimension mismatch! Returned original left hand side." << std::endl;
			return *this;
		}
	}
// =================================================================================================
#pragma endregion

#pragma region GETTER_SETTER
// ========================================DOKUMENTATION============================================
///\fn		ref_fed
///\brief	Gibt die Adresse des ersten Elements in '_data' zurueck.
// -------------------------------------------------------------------------------------------------
	dtype* ref_fed() { return &(_data[0]); }
// =================================================================================================

// ========================================DOKUMENTATION============================================
///\fn		ref_fec
///\brief	Gibt die Adresse des ersten Elements in '_column' zurueck.
// -------------------------------------------------------------------------------------------------
	int* ref_fec() { return &(_col_idx[0]); }
// =================================================================================================

// ========================================DOKUMENTATION============================================
///\fn		ref_fer
///\brief	Gibt die Adresse des ersten Elements in '_row' zurueck.
// -------------------------------------------------------------------------------------------------
	int* ref_fer() { return &(_row_idx[0]); }
// =================================================================================================
#pragma endregion

#pragma region MEMBER_FUNCTIONS
// ========================================DOKUMENTATION============================================
///\fn		printf
///\brief	Ausgabe der Daten als Matrix in der Kommandozeile.
// -------------------------------------------------------------------------------------------------
	void printf()
	{
		for (int i = 0; i < _dim; i++)
		{
			for (int j = 0; j < _dim; j++) std::cout << (*this)(j, i) << " ";
			std::cout << std::endl;
		}
	}
// =================================================================================================
#pragma endregion
   };

