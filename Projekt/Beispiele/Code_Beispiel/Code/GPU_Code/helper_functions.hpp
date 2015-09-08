// ################################################################################################
//				  Header: helper_functions
// ------------------------------------Doxygen-Dokumentation---------------------------------------
///  \file helper_functions.hpp
///  \brief
///  unterstützende Funktionen bzgl. Matrix Operationen
// ------------------------------------------------------------------------------------------------
// ################################################################################################

#include<iostream>
using namespace std;

// ========================================DOKUMENTATION============================================
//					      random_float
///\brief	liefert zufällige Zahl aus (0,1)
// -------------------------------------------------------------------------------------------------
// =================================================================================================   
float random_float()
{
    return rand() / (RAND_MAX + 1.);
}


// ========================================DOKUMENTATION============================================
//					      solutuion_check
///\brief	überprüft das Ergebniss anhand einer zufälligen Zeile TODOOOOO
// -------------------------------------------------------------------------------------------------
///@param	matrix    Matrix
///@param	vector    Vektor
///@param solution  Lösung  
// ================================================================================================= 
bool solution_check()
{
     return true;
}




// ========================================DOKUMENTATION============================================
//					      dense_vec_fill
///\brief	füllt jeden Eintrag eines Vektors mit einem float aus (0,1)
// -------------------------------------------------------------------------------------------------
///@param	vector    1D Array der Größe dim, wird während der Laufzeit gefüllt
///@param	n	     Dimension des Vektor

// =================================================================================================    
void dense_vec_fill(float* vector,int n)
{
     for(int k=0;k<n;k++)
     {
          vector[k]=random_float();
     }

}

// ========================================DOKUMENTATION============================================
//					      sparse_mat_fill
///\brief	füllt ungefähr "1/prozent" zufällige Einträge der Matrix mit einem float Wert zwischen (0,1)
// -------------------------------------------------------------------------------------------------
///@param	matrix    1D Array der Größe dim*dim, wird während der Laufzeit gefüllt
///@param prozent   gibt an welche Anteile der Matrix gefüllt werden sollen
///@param	n	     Dimension der Matrix
// =================================================================================================     
void sparse_mat_fill(float* matrix,int prozent, int n)
{
	for(int i=0;i<n*n;i++)
	{
		int small=rand()%prozent;
		if(small==0)
		{
	          matrix[i]=random_float();
		}
		else
		{
		     matrix[i]=0.0;
		}
	}
}

// ========================================DOKUMENTATION============================================
//					      print_float
///\brief	gibt ein float Array der Größe dim aus
// -------------------------------------------------------------------------------------------------
///@param	input     1D float Array
///@param	n	     Größe des Arrays
// =================================================================================================  

void print_float(float* input,int n)
{
     for(int i=0;i<n;i++)
     {
          cout << input[i] << " ";
     }
     cout << endl << "-------------------------------------------------------------------" << endl;
}


// ========================================DOKUMENTATION============================================
//					      print_mat
///\brief	gibt ein float Array der Größe dim*dim aus
// -------------------------------------------------------------------------------------------------
///@param	input     1D float Array
///@param	n	     Größe des Arrays
// ================================================================================================= 
void print_mat(float* input,int n)
{
     for(int i=0;i<n;i++)
     {
          for(int j=0;j<n;j++)
          {
               cout << input[j+i*n] << "  " << "\t";
          }
     cout << endl;
     } 
     cout << "-------------------------------------------------------------------" << endl;
     
}


// ========================================DOKUMENTATION============================================
//					      print_int
///\brief	gibt ein int Array der Größe dim aus
// -------------------------------------------------------------------------------------------------
///@param	input     1D int Array
///@param	n	     Größe des Arrays
// ================================================================================================= 
void print_int(int* input,int n)
{
     for(int i=0;i<n;i++)
     {
          cout << input[i] << " ";
     }
     cout << endl << "-------------------------------------------------------------------" << endl;
}




























