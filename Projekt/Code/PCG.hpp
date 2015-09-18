#pragma once

#include "Vector.hpp"
#include "DIA.hpp"
#include <iostream>
#include <math.h>
using namespace std;

// kommt bestimmt bald in Vector.hpp
template<typename data,typename data1, typename data2>
void sp(data& result, const Vector<data1>& v1, const Vector<data2>& v2)
// Anwendungsbsp: Skalarprodukt = sp(v1,v2);
{
	assert( v1.dim()==v2.dim() );

	data a(0);
	for(int i=0;i<v1._dim;i++)
	{
		a += static_cast<data>(v1[i])*static_cast<data>(v2[i]);
	}	
	result = a;
}

template <typename restype, typename mattype, typename vectype>
int CG(Vector<restype>& x, DIA<mattype>& A, Vector<vectype>& b) {
	Vector<restype> r(b.dim());	// alle Eintraege werden bei der Initalisierung 0 gesetzt
	matvec(r, A, x);	// r_0 setzen
	r.skalarmult(-1);
	r.vecadd(b);
		Vector<restype> d(r);	// d_0 = r_0
	Vector<restype> z(b.dim());

		restype rr_old(0);
	restype rr_new(1);
	restype dz(1);
	Vector<restype> xold(b.dim());
	Vector<restype> rold(b.dim());

	cout << fixed << setprecision(3);

	double TOL(0.0000000001);

	int k(0);
	for (; k < b.dim() && rr_new > TOL; ++k) {
		matvec(z, A, d);
		sp(rr_old, r, r);
		sp(dz, d, z);
		xold = x;
		x = d;
		x.skalarmult(rr_old / dz);
		x.vecadd(xold);
		rold = r;
		r = z;
		r.skalarmult(-(rr_old / dz));
		r.vecadd(rold);
		sp(rr_new, r, r);	// should be used for the new alpha1
		d.skalarmult(rr_new / rr_old);
		d.vecadd(r);
	}
	return k;
}

template <typename restype, typename mattype, typename vectype>
int PCG_Jacobi(Vector<restype>& x, DIA<mattype>& A, Vector<vectype>& b) {
	int k(0);
	k = CG(x, A, b);
	return k;
}
