#pragma once
#include <iostream>
#include <math.h>
#include <iomanip>
#include <cassert>
using namespace std;



/////////// Ausgelagerte Vektorfunktionen ///////////

// Ausgabe von Vektoren, Ueberladung des "<<"-Operators
template<typename data>
ostream& operator<< (ostream& out, const Vector<data>& v)
// Anwendungsbsp: cout << v << endl;
{
	for(int i=0;i<v._dim;i++)
	{
		cout<< setw(3);	// Abstand zwischen Eintraegen
		out<<v[i]<<" ";
	}
	return out;
};

// Skalarprodukt von zwei Vektoren v1 und v2 (mit Konvertierung des Datenformates (static_cast))
template<typename data,typename data1, typename data2>
data sp(const Vector<data1>& v1, const Vector<data2>& v2)
// Anwendungsbsp: Skalarprodukt = sp(v1,v2);
{
	assert( v1.dim()==v2.dim() );

	data a(0);
	for(int i=0;i<v1._dim;i++)
	{
		a += static_cast<data>(v1[i])*static_cast<data>(v2[i]);
	}
	return a;
}

//Standard 2-Norm: Betrag = norm(x)
template<typename data>
double norm(Vector<data>& y)
{
	double norm2=0;
	for (int i=0;i<y._dim; i++)
	{
		norm2+=y[i]*y[i];
	}
	return sqrt(norm2);							//Ziehe Wurzel, gebe Norm aus
}
#include "Vector.tpp"
