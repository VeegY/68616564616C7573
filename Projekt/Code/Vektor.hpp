#pragma once
#include <iostream>
#include <math.h>
#include <iomanip>
#include <cassert>
using namespace std;

// Implementierung der Vektorklasse auf Basis von Arrays
template<typename data>
class Vector
{
public:
	// Variablen, auf die man global zugreifen kann
	int _dim;	// Dimension des Vektors
	data* _daten;	// Eintraege des Vektors

	// Standardkonstruktor
	Vector(): _daten(NULL)
	// Anwendungsbsp: Vector<double> Vektor1();
	{
		#ifdef CTOR_Test
			std::cout << "Standard-CTOR, no memory allocation " << std::endl;
		#endif		
		
		_dim=0;
	}

	// Konstruktor eines Vektors bestimmter Laenge 
	Vector(const int dim) : _dim(dim)
	// Anwendungsbsp: Vector<double> Vektor2(5);
	{
		#ifdef CTOR_Test
			std::cout << "0-Vektor der Laenge: " << dim << " erstellen" <<  std::endl;
		#endif
        
		_daten = new data[_dim];
		for(int i=0;i<_dim;i++)
		{
				_daten[i]=data(0);
		}
    	}

	// Copy-Konstruktor
	Vector(const Vector& v)
	// Anwendungsbsp: Vector<double> Vektor3(v);
	{
		#ifdef CTOR_Test
			std::cout << "copy-CTOR, memory allocation" << std::endl;
		#endif

		_dim=v.dim();
		_daten = new data[_dim];
		for (int i=0; i < _dim; ++i) {
				_daten[i] = v._daten[i];
		}
	}

	// Copy-Konstruktor mit konvertierung des Datentyps
	template<typename data2>
	Vector(const Vector<data2>& v) 
	// Anwendungsbsp: Vector<float> Neuer_vektor(v);
	{
		#ifdef CTOR_Test
			std::cout << "copy-CTOR (format conversion), memory allocation" << std::endl;
		#endif		

		_dim=v.dim();
		_daten = new data[_dim];
		for (int i=0; i < _dim; ++i) {
				_daten[i] = static_cast<data>(v._daten[i]);
		}
	}

	// Destruktor um den Speicher freizugeben
	~Vector()
	{
		#ifdef CTOR_Test
			std::cout << "DTOR" << std::endl;
		#endif

		delete[] _daten; //Destruktor funktioniert ja nur wenn ich vorher beim CTOR mit new Speicher fuer die Daten freigemacht habe
	}

	// Zugriff auf die Dimension eines Vektors
	int dim() const
	// Anwendungsbsp: dimension = Vector.dim();
	{
		return _dim;
	}

	// Ueberladung des "[]"-Operators (konstant)
    	const data& operator[] (const int index) const 
	// Anwendungsbsp: a = Vector[i]; gibt den i-ten Eintrag des Vektors wieder uns speichert ihn in a
	{
        	return _daten[index];
	}
	
	// Ueberladung des "[]"-Operators (veraenderlich)
 	data& operator[] (const int index) 
	// Anwendungsbsp: Vector[i] = b; veraendert den i-ten Eintrag des Vektors und setzt ihn als b
	{
        	return _daten[index];
	}
	
	// Ueberladung des "="-Operators
	Vector& operator= (const Vector&  v)
	// Anwendungsbsp: Vector Neuervektor = Altervektor;
	{
		assert( _dim == v.dim() );

		for(int i=0;i<_dim;i++)
		{
			_daten[i]=v._daten[i];
		}
		return *this;  
	}

	// Multiplikation mit einem Skalar (mit Konvertierung des Datenformates (static_cast))
	template<typename data2>
	Vector& skalarmult (const data2& skalar)
	// Anwendungsbsp: Vector.skalarmult(skalar);
	{
		for(int i=0;i<_dim;i++)
		{
			_daten[i]*=static_cast<data>(skalar);
		}
		return *this;
	}
	
	// Vektoraddition (gleiche Dimension!) (mit Konvertierung des Datenformates (static_cast))
	template<typename data2>
	Vector& vekadd (const Vector<data2>& v)
	// Anwendungsbsp: Vector.vekadd(v);
	 {
		assert( _dim==v.dim() );

		for(int i=0;i<_dim;i++)
		{
			_daten[i]+=static_cast<data>(v._daten[i]);
		}
		return *this;
	 }

	// Vektorsubtraktion (gleiche Dimension!), Ueberladung des "-="-Operators (mit Konvertierung des Datenformates (static_cast))
	template<typename data2>
	Vector& veksub (const Vector<data2>& v)
	// Anwendungsbsp: Vector.veksub(v);
	{
		assert( _dim==v.dim() );

		for(int i=0;i<_dim;i++)
		{
			_daten[i]-=static_cast<data>(v._daten[i]);
		}
		return *this;
	}
};

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
