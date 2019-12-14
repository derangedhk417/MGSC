#include <iostream>
#include <mathimf.h>
#include <chrono>
#include "gaussian.cpp"
#include <gsl/gsl_min.h>

using namespace std;

int main(int argc, char ** argv) {

	float min = 1.0;
	float max = 1e7;
	int   N   = 10;

	float increment = (max - min) / N;

	cout << "Scale, Error, Expectation" << endl;	

	for (int i = 0; i < N; ++i) {
		// Single Term Testing
		float scale = min + increment * i;
		float A1[]  = {1.0 * scale, 1.0 * scale, 1.0 * scale};
		float ** A  = new float*[1];
		A[0]        = A1;

		float s1[] = {0.0, 0.0, 0.0};
		float ** s = new float*[1];
		s[0]       = s1;

		float C[]  = {1.0};
		float R1[] = {0.0, 0.0, 0.0};
		float ** R = new float*[1];
		R[0]       = R1;

		float Q[] = {N_qe};
		
		GaussianWavefunction wavefn(1, 0.01, 256);
		wavefn.A  = A;
		wavefn.s  = s;
		wavefn.C  = C;
		wavefn.R  = R;
		wavefn.Nu = 1;
		wavefn.Q  = Q; 

		float err;

		auto begin        = chrono::high_resolution_clock::now();
		float expectation = wavefn.getHamiltonianExpectation(&err);
		auto end          = chrono::high_resolution_clock::now();
		auto duration     = chrono::duration_cast<chrono::nanoseconds>(end - begin).count();

		cout << scale << ", " << err << ", " << expectation << endl;
	}


	

	return 0;
}