// This file defines the full gaussian wavefunction and 
// all code necessary to work with it.
#include <iostream>
#include <mathimf.h>
#include <chrono>
#include "constants.h"
#include "integrate.cpp"

// This is the width around the center of the gaussian terms
// for which the integral will be carried out. This is based
// on the assumption that the integral doesn't need to be 
// carried out to infinity, because the integrand will go to
// zero very quickly at distances far from the center of the
// gaussian and the nucleus.
const float range = 5.0;

class GaussianWavefunction {
public:
	int      m;  // Number of terms
	int      n;  // Number of inputs
	float ** A;  // Shape matrices (really vectors due to simplifications)
	float ** s;  // Shift matrices (really vectors)
	float  * C;  // Constants on each Gaussian term
	float ** R;  // List of nuclear coordinates
	int      Nu; // Number of nuclei
	float  * Q;  // List of nuclear charges
	
	float * lowerBounds;
	float * upperBounds;


	DynamicIntegrator integrator;

	GaussianWavefunction(int nterms) {
		// In a future version, this will be 3 times and number of electrons.
		// Right now it is just going to be 3, corresponding to a single 
		// electron.
		n = 3;
		m = nterms;

		
		lowerBounds = new float[3];
		upperBounds = new float[3];

		// This is used to numerically integrate the electon - nucleus
		// integral. We need to reset the bounds each time we use
		// the integrator, so we will just pass the uninitialized
		// arrays for now.
		integrator = new DynamicIntegrator(
			integrandWrapper, lowerBounds, upperBounds, 3
		);
	}

// Math
public:
	// This calculates the normalization constant for the entire 
	// wavefunction.
	float getNormalizationConstant() {
		// We basically need to sum the integral over all space of each
		// term. Each term is multiplied by it's corresponding C coefficient.
		// We then return the reciprocal of this sum.

		float totalsum = 0.0;
		for (int outer_tidx = 0; outer_tidx < m; ++outer_tidx) {
			for (int tidx = 0; tidx < m; ++tidx) {
				float termsum = 1.0;
				for (int iidx = 0; iidx < n; ++iidx) {
					float _s      = s[tidx][iidx] + s[outer_tidx][iidx];
					float _A      = A[tidx][iidx] + A[outer_tidx][iidx];
					float ssquare = _s*_s;
					float aterm   = _A
					termsum *= sqrtf(pi / aterm)*expf(ssquare / (4.0 * aterm));
				}
				termsum *= C[tidx]*C[outer_tidx];
				totalsum += termsum;
			}
		}		

		return 1.0 / totalsum;
	}

	// Gets the expectation value of the Hamiltonian under the current
	// conditions.
	float getHamiltonianExpectation() {
		float B = getNormalizationConstant();
		float T = getKineticExpectation(B);
		float V = getPotentialExpectation(B);

		return T + V;
	}

	// Given the appropriate normalization constant, returns the
	// expectation value of the kinetic energy for the given set 
	// of parameters.
	float getKineticExpectation(float B) {
		float coefficient = -((B*B)*(N_hbar*N_hbar) / (2*N_me))*powrf(pi, (n + 1) / 2.0);
		float sum         = 0.0;

		for (int j_prime = 0; j_prime < m; ++j_prime) {
			for (int j = 0; j < m; ++j) {
				for (int i = 0; i < n; ++i) {
					float D_ijjp       = A[j_prime][i] + A[j][i];
					float F_ijjp       = s[j_prime][i] + s[j][i];
					float coefficient1 = C[j_prime] * C[j];

					float Gamma_t2  = (s[j][i]*s[j][i] - 2*A[j][i])*(1/sqrtf(D_ijjp))
					Gamma_t2       += (A[j][i]*A[j][i])*(2*D_ijjp + F_ijjp*F_ijjp)*(1 / powrf(D_ijjp, 2.5))
					Gamma_t2       -= 2*A[j][i]*s[j][i]*F_ijjp*(1 / powrf(D_ijjp, 1.5))
					float Gamma     = expf((F_ijjp*F_ijjp) / (4*D_ijjp))*Gamma_t2;

					// Now we compute the product term.
					float product = 1.0;
					for (int l = 0; l < n; ++l) {
						if (l != i) {
							float D_ljjp       = A[j_prime][l] + A[j][l];
							float F_ljjp       = s[j_prime][l] + s[j][l];
							product *= (1 / sqrtf(D_ljjp))*expf((F_ljjp*F_ljjp) / (4 * D_ljjp));
						}
					}

					sum += coefficient1*Gamma*product;
				}
			}
		}

		return coefficient * sum;
	}

	// Given the appropriate normalization constant, returns the
	// expectation value of the potential energy for the given set 
	// of parameters.
	float getPotentialExpectation(float B) {
		
	}

// Helper Functions
public:

	// This uses the center location of the given term to set the
	// bounds for the integration intelligently. The terms are 
	// zero indexed, so idx = 0  would correspond to the first
	// term.
	// 
	// tidx = index of gaussian term
	// ridx = index of nucleus
	void initializeBounds(int tidx, int ridx) {
		lowerBounds[0] = fmin(R[ridx][0], -s[tidx][0]) - range;
	    upperBounds[0] = fmax(R[ridx][0], -s[tidx][0]) + range;

	    lowerBounds[1] = fmin(R[ridx][1], -s[tidx][1]) - range;
	    upperBounds[1] = fmax(R[ridx][1], -s[tidx][1]) + range;

	    lowerBounds[1] = fmin(R[ridx][2], -s[tidx][2]) - range;
	    upperBounds[1] = fmax(R[ridx][2], -s[tidx][2]) + range;
	}

	static float integrandWrapper(float * input, void * self) {

	}

};