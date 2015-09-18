// #################################################################################################
//					 FD - Code Projekt
// #################################################################################################
// 				     CPP-File: vector_CPU_TEST
// ------------------------------------Doxygen-Dokumentation----------------------------------------
///  \file vector_CPU_TEST.cpp
///  \brief
///  Testet die Funktionalitaet der Implementierungen in vector_CPU.hpp.
// -------------------------------------------------------------------------------------------------
// Kompilieren mit dem folgenden Compilerbefehl: !! WIP !!
//   gcc vector_CPU_TEST.cpp
// -------------------------------------------------------------------------------------------------
// Verwendete Header:
#	include <iostream>

#	include "vector_metrics.hpp"
// -------------------------------------------------------------------------------------------------
// Copyright : 2015 malte-christian.loddoch@tu-dortmund.de
// #################################################################################################

// =================================================================================================
//					 main-Funktion
// -------------------------------------------------------------------------------------------------
///  \brief Fuehrt den Test der vector_CPU Implementierung durch.
// =================================================================================================
int main(int argc, char** argv)
{
// =================================================================================================
// Initialisiere TEST-Variablen und setze Testbedingungen auf.
// =================================================================================================
	vector<int> test1(3);
	vector<int> test2(3);
	vector<int> test3(3);

	vector<int> coeffs(2);

	test1[0] = 1;
	test1[1] = 0;
	test1[2] = 0;
	test2[0] = 0;
	test2[1] = 1;
	test2[2] = 2;

	coeffs[0] = 3;
	coeffs[1] = 5;

	int* testbasis[2] = { test1.ref_fed(), test2.ref_fed() };
	
// =================================================================================================
// TEST - Copy-Konstruktor
// =================================================================================================
	std::cout << "TEST - copy-constructor:" << std::endl;
	std::cout << "original vector:" << std::endl;
	test1.printf();
	vector<int> cpy_test1(test1);
	std::cout << "copied vector:" << std::endl;
	cpy_test1.printf();
	std::cout << std::endl;

// =================================================================================================
// TEST - '[]'-Operator
// =================================================================================================
	std::cout << "TEST - '[]'- operator:" << std::endl;
	std::cout << "test1:" << std::endl;
	test1.printf();
	std::cout << "test2:" << std::endl;
	test2.printf();
	std::cout << "test1-values: [0], [1], test2-values: [1]" << std::endl;
	std::cout << test1[0] << " " << test1[1] << " " << test2[1] << std::endl << std::endl;

// =================================================================================================
// TEST - Assignment-Operator
// =================================================================================================
	std::cout << "TEST - assignment-operator:" << std::endl;
	std::cout << "cpy_test1 before:" << std::endl;
	cpy_test1.printf();
	cpy_test1 = test2;
	std::cout << "cpy_test1 after:" << std::endl;
	cpy_test1.printf();
	std::cout << std::endl;

// =================================================================================================
// TEST - 'linco'-Funktion
// =================================================================================================
	std::cout << "TEST - linco function:" << std::endl;
	test3.linco(&coeffs, testbasis, 2);
	std::cout << "coefficients:" << std::endl;
	coeffs.printf();
	std::cout << "test1:" << std::endl;
	test1.printf();
	std::cout << "test2:" << std::endl;
	test2.printf();
	std::cout << "result:" << std::endl;
	test3.printf();
	std::cout << std::endl;

// =================================================================================================
// TEST - 'dotprod'-Funktion
// =================================================================================================
	std::cout << "TEST - dotprod function:" << std::endl;
	std::cout << "test1:" << std::endl;
	test1.printf();
	std::cout << "test2:" << std::endl;
	test2.printf();
	std::cout << "result: " << dotprod(test1, test2) << std::endl;
	std::cout << std::endl;

// =================================================================================================
// TEST - Ende der Tests
// =================================================================================================

	return 0;
}
