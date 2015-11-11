﻿//=============================================================================================//
//                                      DIA-Matrix Klasse                                      //
//=============================================================================================//
// enthaelt: (int)        _dim        Matrixdimension (quadratisch)                            //
//           (int)        _numDiags   Anzahl der Baender (Band der Dicke 3 sind 3 Baender)     //
//    (Vector<datatype>*) _data       Matrix Eintraege (in der Reihenfolge des _offset Vektors)//
//    (Vector<int>*)      _offsett    Position der Baender in Bezug auf die Hauptdiagonale     //
//=============================================================================================//
// Anmkerkungen:                                                                               //
// - "TODO" kennzeichnet eine Baustelle. Es kann auch nicht gekennzeichnete geben.             //
// - Noch sind alle Variablen oeffentlich. Das wird wahrscheinlich nicht so bleiben.           //
//   Bitte verwandte Funktionen benutzen. (z.B.: benutze 'matrix.dim()' statt 'matrix._dim' )  //
//=============================================================================================//
// Beispiel:  0   1 7 o o      _dim = 4;                                                       //
//              0 o 2 8 o      _numDiags = 3;                                                  //
//                5 o 3 9      _data ~ {0, 0, 5, 6, 1, 2, 3, 4, 7, 8, 9, 0}  (genaue Info      //
//                o 6 o 4 0    _offset ~ {-2, 0, 1}                           vgl Vector.hpp)  //
//=============================================================================================//

#pragma once
#include <iostream>
#include <math.h>
#include <stdexcept>	// fuer: throw exception
#include "vector.hpp"
using namespace std;

template<typename datatype>
class DIA {
public:
//private:
// Variablen
	int _dim;
	int _numDiags;
	Vector<datatype>* _data;
	Vector<int>* _offset;

//public:
// Konstruktoren
	// Standard Konstruktor
	DIA(): _data(), _offset(), _dim(0), _numDiags(0) { }

	// Konstruktor mit Datenuebergabe
	DIA(int dim, int numDiags, Vector<datatype>& data, Vector<int>& offset) {
		_dim = dim;
		_numDiags = numDiags;
		_data = new Vector<datatype>(data);
		_offset = new Vector<int>(offset);
	}

	// Kopier-Konstruktor
	DIA(const DIA& matrix) {
		_dim = matrix._dim;
		_numDiags = matrix._numDiags;
		_data = new Vector<datatype>(*matrix._data);
		_offset = new Vector<int>(*matrix._offset);
	}

	// Kopier-Konstruktor inkl Konvertierung
	// TODO

	// Destruktor
	~DIA() {
		delete _data;
		delete _offset;
	}

// Daten-Ausgaben
	// Gibt die Dimension der Matrix zurueck
	int dim() { return _dim; }

	// Gibt die Anzahl der Baender zurueck
	int numDiags() { return _numDiags; }

	// Gibt einen Zeiger auf den Datenvektor als Array zurueck (Daten lassen sich darueber nicht veraendern)
	// Bsp: const double* matrixbandeintraege = matrix.data();
	const datatype* data() { return _data->_data; }

	// Gibt einen Zeiger auf den offset-Vektor als Array zurueck (Daten lassen sich darueber nicht veraendern)
	// Bsp: const int* bandpositionen = matrix.offset();
	const int* offset() { return _offset->_data; }

	// Gibt den Matrix-Eintrag an (x, y) zurueck
	datatype value(int x, int y) {
		for (int i(0); i < _numDiags; ++i)
			if (y - x == _offset->_data[i])
				return _data->_data[i * _dim + x];
		return datatype(0);
	}

	// Matrixausgabe (maximal bis 10x10)
	void show() {
		if (_dim > 10)
			cout << "I'm so fat. No one wanna see me! :'(" << endl;
		else
			for (int i(0); i < _dim; ++i) {
				for (int j(0); j < _dim; ++j)
					cout << value(i, j) << " ";
				cout << endl;
			}
	}

// Funktionen
	// Ueberprueft, ob da Nullen stehen, wo Nullen stehen muessen
	// Nur zur Not. Endgueltige Funktionen/Algorithmen sollten ohne auskommen!
	bool checkIntact() {
		for (int i(0); i < _numDiags; ++i) {
			if (_offset->_data[i] < 0) {
				for (int j(0); j < -_offset->_data[i]; ++j)
					if (_data->_data[i*_dim+j] != datatype(0))
						return false;
			}
			else if (_offset->_data[i] > 0) {
				for (int j(_dim - _offset->_data[i]); j < _dim; ++j)
					if (_data->_data[i*_dim+j] != datatype(0))
						return false;
			}
		}
		return true;
	}

	// Setzt ohne zu ueberpruefen alle notwendigen Werte = 0
	// Nur zur Not. Endgueltige Funktionen/Algorithmen sollten ohne auskommen!
	void repair() {
		for (int i(0); i < _numDiags; ++i) {
			if (_offset->_data[i] < 0) {
				for (int j(0); j < -_offset->_data[i]; ++j)
					_data->_data[i*_dim+j] = datatype(0);
			}
			else if (_offset->_data[i] > 0) {
				for (int j(_dim - _offset->_data[i]); j < _dim; ++j)
					_data->_data[i*_dim+j] = datatype(0);
			}
		}
	}

// Operatoren
	// Matrix * Vektor ist (als nicht-Memberfkt) ausgelagert (s.u.)
};


// Matrix mal Vektor: Alternativversion. Nur noch numDiags() if- Abfragen (und innere for-Schleife grundsätzlich vektorisierbar).
// Bsp: matvec(result, mat, vec);
template <typename restype, typename mattype, typename vectype>
void matvec(Vector<restype>& result, DIA<mattype>& mat, Vector<vectype>& vec) {
	int start, end;
	if (result._dim != mat._dim || mat._dim != vec._dim) {
		throw invalid_argument(" -Achtung! Dimensionsfehler! (matvec)- ");
	}
	else {
		for (int i(0); i < result.dim(); ++i){
			result[i] = 0;
		}
		for (int j(0); j < mat.numDiags(); ++j) {
			if (mat.offset()[j] <= 0){
				start = -mat.offset()[j];
				end = mat.dim();
			}
			else{
				start = 0;
				end = mat.dim() - mat.offset()[j];
			}
			for (int i(start); i < end; ++i) {
				result[i] += static_cast<restype>((*mat._data)[mat.dim() * j + i]) * static_cast<restype>(vec[i + mat.offset()[j]]);
			}
		}
	}
}


// Defektberechnung
// Bsp: defect(r, mat, rhs, vec);
template <typename datatype>
void defect(Vector<datatype>& r, DIA<datatype>& mat, Vector<datatype>& rhs, Vector<datatype>& vec) {
	int start, end;
	if (r._dim!=mat._dim || mat._dim!=rhs._dim || rhs._dim!=vec._dim) {
		throw invalid_argument(" -Achtung! Dimensionsfehler! (defect)- ");
	}
	else {
		r = rhs;	// da dies hier gesetzt wird, kann man innerhalb der Schleifen "-=" nutzen
		for (int j(0); j < mat.numDiags(); ++j) {
			if (mat.offset()[j] <= 0){
				start = -mat.offset()[j];
				end = mat.dim();
			}
			else{
				start = 0;
				end = mat.dim() - mat.offset()[j];
			}
			for (int i(start); i < end; ++i) {
				r[i] -= (*mat._data)[mat.dim() * j + i] * vec[i + mat.offset()[j]];
			}
		}
	}
}

