#include <iostream>
#include <mathimf.h>
#include <chrono>
#include "gaussian.cpp"

using namespace std;

int main(int argc, char ** argv) {

	float test[] = {
		1e-2, 2e-2, 4e-2,
		1e-1, 2e-2, 4e-1,
		1.0,  2.0,  3.0,
		10.0, 20.0, 30.0,
		40.0, 50.0, 60.0,
		70.0, 80.0, 90.0,
		100.0, 1000.0
	};

	int n_test = 20;

	for (int i = 0; i < n_test; ++i) {
		// Single Term Testing
		float scale = test[i];
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
		
		GaussianWavefunction wavefn(1);
		wavefn.A  = A;
		wavefn.s  = s;
		wavefn.C  = C;
		wavefn.R  = R;
		wavefn.Nu = 1;
		wavefn.Q  = Q; 

		float expectation = wavefn.getHamiltonianExpectation();

		cout << scale << " " << expectation << "eV" << endl;
	}
	

	return 0;
}