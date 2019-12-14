#include <iostream>
#include <mathimf.h>
#include <chrono>
#include "gaussian.h"
#include "constants.h"

using namespace std;

int main() {

	double min = 1.0;
	double max = 1e7;
	int   N    = 10;

	double increment = (max - min) / N;

	cout << "Scale, Expectation" << endl;	

	for (int i = 0; i < N; ++i) {
		// Single Term Testing
		double scale = min + increment * i;
		double A1[]  = {1.0 * scale, 1.0 * scale, 1.0 * scale};
		double ** A  = new double*[1];
		A[0]         = A1;

		double s1[] = {0.0, 0.0, 0.0};
		double ** s = new double*[1];
		s[0]        = s1;

		double C[]  = {1.0};
		double R1[] = {0.0, 0.0, 0.0};
		double ** R = new double*[1];
		R[0]        = R1;

		double Q[] = {N_qe};
		
		GaussianWFN wavefn(1, 1);
		wavefn.A  = A;
		wavefn.s  = s;
		wavefn.C  = C;
		wavefn.R  = R;
		wavefn.Nu = 1;
		wavefn.Q  = Q; 
		
		double expectation = wavefn.getHamiltonianExpectation();
	
		cout << scale << ", " << expectation << endl;
	}
	
	return 0;
}