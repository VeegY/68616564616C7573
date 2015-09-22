#include <iostream>
using namespace std;

#include "DIA.hpp"
#include "Vector.hpp"
#include "PCG.hpp"
#include <time.h>

int main() {
	cout << "this is not a test. this is serious stuff. please stay calm and dont panic." << endl;
	cout << "haha! you idiots! ruun!!" << endl << endl;

	clock_t start;
	double TOL(0.00001);
	double exact(0.0);

	int dim[10] = { 0 }; for(int i(0);i<10;++i) dim[i]=10000;
	dim[3] = 10;
	double time[10] = { 0 };
	int numIt[10] = { 0 };
	int valx[10] = { 0 };
	bool OK[10]; for(int i(0);i<10;++i) OK[i] = true;

	// Test 1 - Einheitsmatrix
	Vector<double> v1(dim[1]); v1.set(1.0);
	Vector<int> o1(1);	// wird 0 gesetzt
	DIA<double> A1(dim[1], o1.dim(), v1, o1);
	Vector<double> b1(dim[1]); b1.set(1.0);
	Vector<double> x1(dim[1]); x1.set(0.0);
	start = clock();
	numIt[1] = CG(x1, A1, b1);
	time[1] = static_cast<double>((clock() - start)) / static_cast<double>(CLOCKS_PER_SEC);
	for (; valx[1] < dim[1] && OK[1]; ++ valx[1])
		if (x1[valx[1]] != 1.0)
			OK[1] = false;
	if (OK[1])
		cout << "Id mat | dim: " << dim[1] << " | iterations: " << numIt[1] << " | time: " << time[1] << "sec" << endl;
	else
		cout << "FAILED test 1 - first incorrect x-value: " << valx[1] << endl;

	// Test 2 - 1D WL-Gl, b=1
	Vector<double> v2(3*dim[2]); v2.set(1.0); for(int i(0);i<dim[2];++i) v2[dim[2]+i]=-2.0;
	Vector<int> o2(3); o2[0] = -1; o2[1] = 0; o2[2] = 1;
	DIA<double> A2(dim[2], o2.dim(), v2, o2);
	Vector<double> b2(dim[2]); b2.set(1.0);
	Vector<double> x2(dim[2]); x2.set(0.0);
	cout << PCG_Jacobi(x2, A2, b2) << endl;	// makes following CG useless/ non-informative
	start = clock();
	numIt[2] = CG(x2, A2, b2);
	time[2] = static_cast<double>((clock() - start)) / static_cast<double>(CLOCKS_PER_SEC);
	exact = 0.0;
	for (; valx[2] < dim[2] && OK[2]; ++valx[2]) {
		if (valx[2] <= dim[2] / 2) exact -= (dim[2] - (2 * valx[2])) / 2;
		else exact += (valx[2] - (dim[2] / 2));
		if (x2[valx[2]] < exact - TOL || x2[valx[2]] > exact + TOL)
			OK[2] = false;
	}
	if (OK[2])
		cout << "1D heat equation, b=1 | dim: " << dim[2] << " | iterations: " << numIt[2] << " | time: " << time[2] << "sec" << endl;
	else
		cout << "FAILED test 2 - first incorrect x-value: " << valx[2] << endl;






















	return 0;
}
