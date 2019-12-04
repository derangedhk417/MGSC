#include <iostream>
#include <mathimf.h>
#include <chrono>
#include "gaussian.cpp"

using namespace std;

int main(int argc, char ** argv) {

	// float errors[] = {
	// 	0.1, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01, 0.005, 0.0025
	// };

	float scales[] = {
		1e-10, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e10
	};

	int n = 13;

	cout << "Scale" << " || " << "Time" << " || " << "Error" << " || " << "Value" << endl;

	for (int i = 0; i < n; ++i) {
		// Single Term Testing
		float scale = scales[i];
		float A1[]  = {1.0 * scale, 1.0 * scale, 1.0 * scale};
		float ** A  = new float*[1];
		A[0]        = A1;

		float s1[] = {0.0, 0.0, 0.0};
		float ** s = new float*[1];
		s[0]       = s1;

		float C[]  = {1.0};
		float R1[] = {10.0, 10.0, 10.0};
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

		cout << scales[i] << " || " << duration / 1e6 << " ms || " << err << " eV || " << expectation << " eV" << endl;
	}
	
	

	return 0;
}