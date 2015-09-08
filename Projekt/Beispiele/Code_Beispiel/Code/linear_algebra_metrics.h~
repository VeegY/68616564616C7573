// ################################################################################################
//					 FD - Code Project
// ################################################################################################
//				  Header: linear_algebra_metrics
// ------------------------------------Doxygen-Dokumentation---------------------------------------
///  \file linear_algebra_metrics.h
///  \brief
///  Hilfsfunktionen aus der linearen Algebra basierend auf der Matrix-Vektor Implementierung.
// ------------------------------------------------------------------------------------------------
// Verwendete Header:
#	pragma once
#	include <cmath>

#	include "vector_CPU.h"
#	include "CSR_CPU.h"
// -------------------------------------------------------------------------------------------------
// Verwendet in :
//   CSR_CPU.h
// -------------------------------------------------------------------------------------------------
// Copyright : 2015 malte-christian.loddoch@tu-dortmund.de
// #################################################################################################

#pragma region VEC_FUNCTIONS
// ========================================DOKUMENTATION============================================
//					      dotprod
///\brief	Wertet das Skalarprodukt zweier Vektoren aus.
// -------------------------------------------------------------------------------------------------
///@tparam	dtype_dotprod 	Template Type fuer die 'dotprod' Funktion.
///@param	a		Erster Vektor im Skalarprodukt.
///@param	b		Zweiter Vektor im Skalarprodukt.
// =================================================================================================
   template<typename dtype_dotprod>
   dtype_dotprod dotprod(vector<dtype_dotprod>& a, vector<dtype_dotprod>& b)
   {
	// get proper pointers to the data, otherwise, compilers cannot vectorise the loop
	const dtype_dotprod* __MYRESTRICT ap = a.ref_fed();
	const dtype_dotprod* __MYRESTRICT bp = b.ref_fed();

	dtype_dotprod result(dtype_dotprod(0.0)); ///< result Ergebnis
	#pragma omp parallel for default(shared) reduction(+:result)
	for (int i = 0; i < a._dim; i++) result += ap[i] * bp[i];

	return result;
   }
// =================================================================================================

// ========================================DOKUMENTATION============================================
//					      norm_p
///\brief	Berechnet die p-Norm eines Vektors.
// -------------------------------------------------------------------------------------------------
///@tparam	dtype 		Template Type fuer die 'norm_p'-Funktion.
///@param	vec_inp		Der auszuwertende Vektor.
///@param	p 		Der Grad der Norm.
///@param	result 		Ergebnis der p-Norm.
// =================================================================================================
   template<typename dtype>
   void norm_p(vector<dtype>& vec_inp, int p, dtype& result)
   {
	const dtype* __MYRESTRICT vec_in = vec_inp.ref_fed();

	result = dtype(0.0);
//	#pragma omp parallel for default(shared) reduction(+:result) -- WIP -- : reduction
	for (int i(0); i < vec_inp._dim; i++) result = result + pow(vec_in[i], p);

	result = pow(result, 1.0 / p);
   }
// =================================================================================================
#pragma endregion

#pragma region CSR_FUNCTIONS
// ========================================DOKUMENTATION============================================
///\fn		template<typename dtype_SpMV> void SpMV(csr<dtype_SpMV>& csr_in, vector<dtype_SpMV>& vec_in, vector<dtype_SpMV>& vec_out)
///\brief	Berechnet das Matrix-Vektor Produkt.
// -------------------------------------------------------------------------------------------------
///@tparam	dtype_SpMV 	Template Type fuer die 'SpMV'-Funktion.
///@param	csr_in 		Eingabematrix
///@param	vec_in		Eingabevektor
///@param	vec_out		Ausgabevektor
// =================================================================================================
   template<typename dtype_SpMV> 
   void SpMV(csr<dtype_SpMV>& csr_in, vector<dtype_SpMV>& vec_in, vector<dtype_SpMV>& vec_out)
   {
	// get proper pointers to the data, otherwise, compilers cannot vectorise the loop
	const dtype_SpMV* __MYRESTRICT data = csr_in.ref_fed();
	const dtype_SpMV* __MYRESTRICT col = csr_in.ref_fec();
	const dtype_SpMV* __MYRESTRICT row = csr_in.ref_fer();
	const dtype_SpMV* __MYRESTRICT x = vec_in.ref_fed();
	dtype_SpMV* __MYRESTRICT y = vec_out.ref_fed();

	// calculate SpMV on the pointers
	#pragma omp parallel for default(shared)
	for (int i = 0; i < csr_in._dim; i++)
	{
		//temporary value
		dtype_SpMV temp(dtype_SpMV(0.0)); ///< temp Temporaere Hilfsvariable
		for (int j = row[i]; j < row[i + 1]; j++)
		{	
			temp += data[j] * x[col[j]];
		}
		y[i] = temp;
	}
   }
// =================================================================================================

// ========================================DOKUMENTATION============================================
///\fn		template<typename dtype_defect> void defect(csr<dtype_defect>& csr_in, vector<dtype_defect>& vec_diff, vector<dtype_defect>& vec_in, vector<dtype_defect>& vec_out, dtype_defect& norm_defect)
///\brief	Berechnet den Defekt
// -------------------------------------------------------------------------------------------------
///@tparam	dtype_defect 	Template Type fuer die 'defect'-Funktion .
///@param	csr_in 		Eingabematrix
///@param	vec_diff	Differenz-Eingabevektor
///@param	vec_in		SpMV-Eingabevektor
///@param	vec_out		Ausgabevektor
///@param	norm_defect	Ergebnis der Defektberechnung
// =================================================================================================
   template<typename dtype_defect>
   void defect(csr<dtype_defect>& csr_in, vector<dtype_defect>& vec_diff,
			vector<dtype_defect>& vec_in, vector<dtype_defect>& vec_out,
			dtype_defect& norm_defect)
   {
	// get proper pointers to the data, otherwise, compilers cannot vectorise the loop
	const dtype_defect* __MYRESTRICT data = csr_in.ref_fed();
	const dtype_defect* __MYRESTRICT col = csr_in.ref_fec();
	const dtype_defect* __MYRESTRICT row = csr_in.ref_fer();
	const dtype_defect* __MYRESTRICT b = vec_diff.ref_fed();
	const dtype_defect* __MYRESTRICT x = vec_in.ref_fed();
	dtype_defect* __MYRESTRICT y = vec_out.ref_fed();

	dtype_defect temp; ///< temp Temporaere Hilfsvariable

	// calculate SpMV on the pointers
	#pragma omp parallel for default(shared) reduction(+:temp)
	for (int i = 0; i<(csr_in._dim); i++)
	{
		temp = dtype_defect(0.0);
		for (int j = (row[i]); j<row[i + 1]; j++)
		{
			temp += (data[j] * static_cast<dtype_defect>(x[col[j]]));
		}
		y[i] = b[i] - temp;
	}
	norm_p(vec_out, 2, norm_defect);
   }
// =================================================================================================
#pragma endregion
