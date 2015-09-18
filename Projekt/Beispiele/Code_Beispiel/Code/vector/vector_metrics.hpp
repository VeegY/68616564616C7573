// ################################################################################################
//					 FD - Code Project
// ################################################################################################
//				  Header: vector_metrics
// ------------------------------------Doxygen-Dokumentation---------------------------------------
///  \file vector_metrics.hpp
///  \brief
///  Hilfsfunktionen aus der linearen Algebra basierend auf der Vektor Implementierung.
// ------------------------------------------------------------------------------------------------
// Verwendete Header:
#	pragma once
#	include <cmath>

#	include "vector_CPU.hpp"
// -------------------------------------------------------------------------------------------------
// Verwendet in :
//   vector_CPU_TEST.hpp
// -------------------------------------------------------------------------------------------------
// Copyright : 2015 malte-christian.loddoch@tu-dortmund.de
// #################################################################################################

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
