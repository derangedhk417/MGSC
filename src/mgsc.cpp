// Author:  Adam J. Robinson
// Summary: This file contains a simple attempt at using a monte-carlo like
//          method to find the ground state energy of the hydrogen atom. It
//          works by simply guessing random terms to add to the wavefunction
//          and then accepting or rejecting them based on whether or not they
//          lower the expectation value.

#include <iostream>
#include <mathimf.h>
#include <chrono>
#include "gaussian.h"
#include "constants.h"
#include "random.h"

using namespace std;

int main() {
	// Initialize a wavefunction with a single term.
	GaussianWFN wavefn(1, 1);

	// Start with the best value of a for a single
	// term. 
	double initial = 3.945310E+06;

	double * A0 = new double[3];
	A0[0] = initial;
	A0[1] = initial;
	A0[2] = initial;

	double ** A = new double*[1];
	A[0]        = A0;
	wavefn.A    = A;

	double * s0 = new double[3];
	s0[0] = 0.0;
	s0[1] = 0.0;
	s0[2] = 0.0;

	double ** s = new double*[1];
	s[0]        = s0;
	wavefn.s    = s;

	double * C = new double[1];
	C[0]       = 1.0;
	wavefn.C   = C;

	double * R0 = new double[3];
	double ** R = new double*[1];
	R[0]        = R0;
	wavefn.R    = R;

	double * Q = new double[1];
	Q[0]       = N_qe;

	wavefn.Nu = 1;
	wavefn.Q  = Q;

	NormalDistribution *A_rng = new NormalDistribution(initial, initial);
	NormalDistribution *C_rng = new NormalDistribution(1.0, 2.0);

	// Stores the last excepted expectation value of the Hamiltonian.
	double lastEXP  = 0.0;
	double lastErr  = 0.0;
	double lastStd  = initial;
	double lastStd2 = 5.0;
	int    nTerms   = 1;   // Current number of accepted terms.
	int    nGuess   = 0;   // Current number of terms guessed.

	lastEXP = wavefn.getHamiltonianExpectation(&lastErr);
	wavefn.setMaxCalls(5000);

	cout << "Guess, Expectation" << endl;
	for (int i = 0; i < 1000; ++i) {
		// Guess a new term to add.
		double * newA  = new double[3];
		double * newS  = new double[3];
		double   newC  = C_rng->read();
		double   A_val = abs(A_rng->read());

		newA[0] = A_val;
		newA[1] = A_val;
		newA[2] = A_val;

		newS[0] = 0.0;
		newS[1] = 0.0;
		newS[2] = 0.0;

		wavefn.pushTerm(newA, newS, newC);

		double exp = wavefn.getHamiltonianExpectation(&lastErr);

		if (exp < lastEXP) {
			cout << nGuess << ", " << exp << endl;
			lastEXP = exp;
			nTerms += 1;
		} else {
			wavefn.popTerm();
		}

		nGuess ++;
	}

	cout << "Nterms: " << nTerms << endl;
	cout << "NGuess: " << nGuess << endl;


	// double min = 1.0;
	// double max = 1e7;
	// int   N    = 512;

	// double increment = (max - min) / N;

	// cout << "Scale, Error, Expectation" << endl;	

	// for (int i = 0; i < N; ++i) {
	// 	// Single Term Testing
	// 	double scale = min + increment * i;
	// 	double A1[]  = {1.0 * scale, 1.0 * scale, 1.0 * scale};
	// 	double ** A  = new double*[1];
	// 	A[0]         = A1;

	// 	double s1[] = {0.0, 0.0, 0.0};
	// 	double ** s = new double*[1];
	// 	s[0]        = s1;

	// 	double C[]  = {1.0};
	// 	double R1[] = {0.0, 0.0, 0.0};
	// 	double ** R = new double*[1];
	// 	R[0]        = R1;

	// 	double Q[] = {N_qe};
		
	// 	GaussianWFN wavefn(1, 1);
	// 	wavefn.A  = A;
	// 	wavefn.s  = s;
	// 	wavefn.C  = C;
	// 	wavefn.R  = R;
	// 	wavefn.Nu = 1;
	// 	wavefn.Q  = Q; 
		
	// 	double error;
	// 	double expectation = wavefn.getHamiltonianExpectation(&error);
	
	// 	cout << scale << ", " << error << ", " << expectation << endl;
	// }
	
	return 0;
}