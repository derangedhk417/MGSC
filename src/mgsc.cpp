#include <iostream>
#include <mathimf.h>
#include <chrono>
#include "gaussian.cpp"

using namespace std;

int main(int argc, char ** argv) {

	// float scales[] = {
	// 	980.0,  990.0, 995.0, 997.5, 999.0,
	// 	1001.0, 1002.5, 1005.0, 1010.0, 1020.0 
	// };

	// int n = 10;

	// for (int i = 0; i < n; ++i) {
		// Single Term Testing
		float scale = 2000.0;
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
		
		GaussianWavefunction wavefn(1, 0.01, 128);
		wavefn.A  = A;
		wavefn.s  = s;
		wavefn.C  = C;
		wavefn.R  = R;
		wavefn.Nu = 1;
		wavefn.Q  = Q; 

		float expectation = wavefn.getHamiltonianExpectation();

		cout << scale << " " << expectation << " eV" << endl;
	//}
	
	

	return 0;
}