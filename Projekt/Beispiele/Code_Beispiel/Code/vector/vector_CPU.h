// #################################################################################################
//					 FD - Code Projekt
// #################################################################################################
// 					 Header: vector_CPU
// ------------------------------------Doxygen-Dokumentation----------------------------------------
///  \file vector_CPU.h
///  \brief
///  CPU Implementierung der Vektorklasse unter Verwendung von 'by Reference' Adressierung
///  fuer die  algebraischen Funktionen.
// -------------------------------------------------------------------------------------------------
// Verwendete Header:
#    pragma once
#    include <omp.h>

#    include "Miscellaneous.h"
// -------------------------------------------------------------------------------------------------
// Verwendet in :
//   vector_CPU_TEST.cpp
//   CSR_CPU.h
// -------------------------------------------------------------------------------------------------
// Copyright : 2015 malte-christian.loddoch@tu-dortmund.de
// #################################################################################################

// ========================================DOKUMENTATION============================================
///\class 	vector
///\brief 	Implementiert die Vektorklasse.
// -------------------------------------------------------------------------------------------------
///@tparam	dtype 	Template Type fuer die 'vector'-Klasse
// =================================================================================================
   template<typename dtype> ///
   class vector
   {
public:

#pragma region PUBLIC_VARIABLES
	int _dim;	///<	_dim 	Dimension des Vektors.
	dtype* _data;	///<	_data	Array der Vektoreintraege.
#pragma endregion

#pragma region CONSTRUCTORS
// ========================================DOKUMENTATION============================================
//				       Standard-Konstruktor
///\brief	Definition des Standard-Konstruktor.
// -------------------------------------------------------------------------------------------------
	vector() : _data(NULL) { _dim = 0; }
// =================================================================================================

// ========================================DOKUMENTATION============================================
//					  Dim-Konstruktor
///\brief	Definition eines Konstruktors fuer waehlbare Vektorlaengen.
// -------------------------------------------------------------------------------------------------
///@param 	dim		Groesse des Array '_data'.
// -------------------------------------------------------------------------------------------------
	vector(int dim) : _dim(dim)
	{
		_data = new dtype[dim];
		#pragma omp parallel for default(shared)
		for (int i = 0; i < dim; i++) { _data[i] = dtype(0); }

	}
// =================================================================================================

// ========================================DOKUMENTATION============================================
//					 Copy-Konstruktor
///\brief	Definition des Copy Konstructors
// -------------------------------------------------------------------------------------------------
///@param 	init_vector	Addresse des zu kopierenden 'vector'-Objekts.
// -------------------------------------------------------------------------------------------------
	vector(vector<dtype>& init_vector) : _dim(init_vector._dim)
	{   
		const dtype* __MYRESTRICT data = init_vector.ref_fed();

		_data = new dtype[_dim];
		#pragma omp parallel for default(shared)
		for (int i = 0; i < _dim; i++) { _data[i] = data[i]; }
	}
// =================================================================================================

// ========================================DOKUMENTATION============================================
//					    Destruktor
///\brief	Definition des Destruktors
// -------------------------------------------------------------------------------------------------
	~vector() { delete[] _data; }
// =================================================================================================
#pragma endregion

#pragma region OPERATORS
// ========================================DOKUMENTATION============================================
//					   '[]'-Operator
///\brief	Ueberladung des '[]'-Operator
// -------------------------------------------------------------------------------------------------
///@param	index		Position dessen zu addressierenden Elements.
// -------------------------------------------------------------------------------------------------
	dtype& operator[] (const int index) { return _data[index]; }
// =================================================================================================

// ========================================DOKUMENTATION============================================
//					Assignment-Operator
///\brief	Ueberladung des Assignment-Operator
// -------------------------------------------------------------------------------------------------
///@param	vec		Referenz des Vektors dessen Daten zugewiesen werden sollen.
// -------------------------------------------------------------------------------------------------
	vector<dtype>& operator= (const vector<dtype>& vec)
	{
		// check if Dimensions fit
		if (_dim == vec._dim)
		{
			// allocate memory
			_data = new dtype[_dim];

			// match values
			#pragma omp parallel for default(shared)
			for (int i = 0; i < _dim; i++) _data[i] = vec._data[i];

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
// ========================================DOKUMENTATION===========================================
///\fn		ref_fed
///\brief	Gibt die Adresse des ersten Elements in _data zurueck.
// -------------------------------------------------------------------------------------------------
	dtype* ref_fed() { return &(_data[0]); }
// =================================================================================================
#pragma endregion

#pragma region MEMBER_FUNCTIONS
// ========================================DOKUMENTATION============================================
///\fn		printf
///\brief	Ausgabe der Daten als Zeilenvektor in der Kommandozeile.
// =================================================================================================
	void printf()
	{
		for (int i = 0; i < _dim; i++) std::cout << _data[i] << " ";
		std::cout << std::endl;
	}
// =================================================================================================

// ========================================DOKUMENTATION============================================
///\fn		template<typename dtype_set_all> void set_all(dtype_set_all value)
///\brief	Setzt alle Werte des Vektors auf den Wert 'value'.
// -------------------------------------------------------------------------------------------------
///@tparam	dtype_set_all 	Template Type fuer die 'set_all' Funktion.
///@param	value		Der zu setzende Wert.
// =================================================================================================
	template<typename dtype_set_all> 
	void set_all(dtype_set_all value)
	{
		#pragma omp parallel for default(shared)
		for (int i = 0; i < _dim; i++) _data[i] = value;
	}
// =================================================================================================

// ========================================DOKUMENTATION============================================
///\fn		template<typename dtype_linco> void linco(vector<dtype_linco>* coeffs, dtype_linco** basis, int count)
///\brief	Wertet eine Linearkombination verschiedener Vektoren aus.
// -------------------------------------------------------------------------------------------------
///@tparam	dtype_linco 	Template Type fuer die 'linco' Funktion.
///@param	coeffs		Referenz des Koeffizienten vector.
///@param	basis		Ein Array mit Referenzen zu Vektordaten.
///@param	count		Anzahl Vektor in 'basis'.
// =================================================================================================
	template<typename dtype_linco>
	void linco(vector<dtype_linco>* coeffs, dtype_linco** basis, int count)
	{
		// get proper pointers to the data, otherwise, compilers cannot vectorise the loop
		const dtype_linco* __MYRESTRICT coeffsp = coeffs->ref_fed();

		for (int j = 0; j < count; j++)
		{
			// get proper pointers to the data, otherwise, compilers cannot vectorise the loop
			const dtype_linco* __MYRESTRICT basisp = basis[j];

			#pragma omp parallel for default(shared)
			for (int i = 0; i < _dim; i++) _data[i] += coeffsp[j] * basisp[i];
		}
	}
// =================================================================================================
#pragma endregion
   };
