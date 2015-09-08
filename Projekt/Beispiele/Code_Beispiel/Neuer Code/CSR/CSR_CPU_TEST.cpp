// #################################################################################################
//					 FD - Code Projekt
// #################################################################################################
// 				     CPP-File: CSR_CPU_TEST
// ------------------------------------Doxygen-Dokumentation----------------------------------------
///  \file CSR_CPU_TEST.cpp
///  \brief
///  Testet die Funktionalitaet der Implementierungen in CSR_CPU.h.
// -------------------------------------------------------------------------------------------------
// Kompilieren mit dem folgenden Compilerbefehl: !! WIP !!
//   gcc CSR_CPU_TEST.cpp
// -------------------------------------------------------------------------------------------------
// Verwendete Header:
#	include <iostream>

#	include "linear_algebra_metrics.h"
// -------------------------------------------------------------------------------------------------
// Copyright : 2015 malte-christian.loddoch@tu-dortmund.de
// #################################################################################################

// =================================================================================================
//  					  main-Funktion
// -------------------------------------------------------------------------------------------------
///  \brief Fuehrt den Test der CSR_CPU Implementierung durch. </summary>
// =================================================================================================*/
int main(int argc, char** argv)
{
// =================================================================================================
// Initialisiere TEST-Variablen und setze Testbedingungen auf.
// =================================================================================================
	vector<int> test1(3);
	vector<int> test2(3);
	vector<int> test3(3);

	test1[0] = 1;
	test1[1] = 0;
	test1[2] = 0;
	test2[0] = 0;
	test2[1] = 1;
	test2[2] = 0;

	int norm_defect(0);
	int data1[5] = { 1, 2, 1, 1, 1 };
	int data2[5] = { 3, 1, 3, 6, 1 };
	int row1[4] = { 0, 2, 4, 5 };
	int row2[4] = { 0, 1, 4, 5 };
	int col1[5] = { 0, 1, 0, 1, 2 };
	int col2[5] = { 0, 0, 1, 2, 2 };
	csr<int> mat1(3, 5, data1, col1, row1);
	csr<int> mat2(3, 5, data2, col2, row2);

// =================================================================================================
// TEST - copy-constructor
// =================================================================================================
	std::cout << "TEST - copy-constructor:" << std::endl;
	std::cout << "original matrix:" << std::endl;
	mat1.printf();
	csr<int> cpy_matrix(mat1);
	std::cout << "copied matrix:" << std::endl;
	cpy_matrix.printf();
	std::cout << std::endl;

// =================================================================================================
// TEST - '()'-Operator 
// =================================================================================================
	std::cout << "TEST - '()'- operator:" << std::endl; 
	std::cout << "values: (1,0), (0,2), (1,1)" << std::endl;
	std::cout << mat1(1, 0) << " " << mat1(0, 2) << " " << mat1(1, 1) << std::endl << std::endl;
	
// =================================================================================================
// TEST - assg-Operator
// =================================================================================================
	std::cout << "TEST - assignment-operator:" << std::endl;
	std::cout << "mat2 before:" << std::endl;
	mat2.printf();
	mat2 = mat1;
	std::cout << "mat2 after:" << std::endl;
	mat2.printf();
	std::cout << std::endl;

// =================================================================================================
// TEST - printf-Funktion
// =================================================================================================
	std::cout << "TEST - printf function:" << std::endl;
	std::cout << "matrix:" << std::endl;
	mat1.printf();
	std::cout << std::endl;

// =================================================================================================
// TEST - SpMV-Funktion
// =================================================================================================
	std::cout << "TEST - SpMV function:" << std::endl;
	SpMV(mat1, test2, test3);
	std::cout << "matrix:" << std::endl;
	mat1.printf();
	std::cout << "test1:" << std::endl;
	test1.printf(); 
	std::cout << "result:" << std::endl;
	test3.printf();
	std::cout << std::endl;

// =================================================================================================
// TEST - Defekt-Funktion
// =================================================================================================
	std::cout << "TEST - defect function:" << std::endl;
	defect(mat1, test1, test2, test3, norm_defect);
	std::cout << "matrix:" << std::endl;
	mat1.printf();
	std::cout << "test1:" << std::endl;
	test1.printf();
	std::cout << "test2:" << std::endl;
	test2.printf();
	std::cout << "result:" << std::endl;
	test3.printf();
	std::cout << std::endl;

// =================================================================================================
// TEST - end
// =================================================================================================
	std::cin.get();

	return 0;
}
