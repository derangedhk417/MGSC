// This file defines the full gaussian wavefunction and 
// all code necessary to work with it.
#include <iostream>
#include <mathimf.h>
#include <chrono>
#include "constants.h"
//#include "en_integrator.cpp"
//#include "en_integrator_alt.cpp"
#include "en_integrator_monte_carlo.cpp"

#define DEBUG
#undef  DEBUG

// This is the width around the center of the gaussian terms
// for which the integral will be carried out. This is based
// on the assumption that the integral doesn't need to be 
// carried out to infinity, because the integrand will go to
// zero very quickly at distances far from the center of the
// gaussian and the nucleus.
const float range = 10.0;

// This is passed to the integrand for the electron - nucleus
// integral by the DynamicIntegrator class.
struct IntegrandParameters {
	float A1;
	float A2;
	float A3;
	float s1;
	float s2;
	float s3;
	float *R;
};

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

	double relerr;
	int    size;

	GaussianWavefunction(int nterms, double relerr, int size) {
		// In a future version, this will be 3 times and number of electrons.
		// Right now it is just going to be 3, corresponding to a single 
		// electron.
		n = 3;
		m = nterms;

		this->relerr = relerr;
		this->size   = size;
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
					float aterm   = _A;
					termsum *= sqrtf(pi / aterm)*expf(ssquare / (4.0 * aterm));
				}
				termsum *= C[tidx]*C[outer_tidx];
				totalsum += termsum;
			}
		}		

		return 1.0 / sqrtf(totalsum);
	}

	// Gets the expectation value of the Hamiltonian under the current
	// conditions.
	float getHamiltonianExpectation(float *err) {
		float B = getNormalizationConstant();
		float T = getKineticExpectation(B);
		float V = getPotentialExpectation(B, err);


#ifdef DEBUG

		std::cout << "B: " << B << std::endl;
		std::cout << "T: " << T << std::endl;
		std::cout << "V: " << V << std::endl;

#endif

		return T + V;
	}

	// Given the appropriate normalization constant, returns the
	// expectation value of the kinetic energy for the given set 
	// of parameters.
	float getKineticExpectation(float B) {
		float coefficient = -((B*B)*(N_hbar*N_hbar) / (2*N_me))*powrf(pi, ((float)n) / 2.0);
		float sum         = 0.0;

		for (int j_prime = 0; j_prime < m; ++j_prime) {
			for (int j = 0; j < m; ++j) {
				for (int i = 0; i < n; ++i) {
					float D_ijjp       = A[j_prime][i] + A[j][i];
					float F_ijjp       = s[j_prime][i] + s[j][i];
					float coefficient1 = C[j_prime] * C[j];

					float Gamma_t2  = (s[j][i]*s[j][i] - 2*A[j][i])*(1/sqrtf(D_ijjp));
					Gamma_t2       += (A[j][i]*A[j][i])*(2*D_ijjp + F_ijjp*F_ijjp)*(1 / powrf(D_ijjp, 2.5));
					Gamma_t2       -= 2*A[j][i]*s[j][i]*F_ijjp*(1 / powrf(D_ijjp, 1.5));
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
	float getPotentialExpectation(float B, float * err) {
		float coefficient = (N_qe*B*B) / N_fpien;

#ifdef DEBUG
		std::cout << "Outer Coefficient: " << coefficient << std::endl;
#endif

		float sum = 0.0;
		for (int iu = 0; iu < Nu; ++iu) {
			// Normally there would be another loop here over 
			// the number of electrons, but I'm omitting that in this
			// version because of time constraints. I am only doing
			// one electron systems for now.

			for (int w = 0; w < m; ++w) {
				for (int l = 0; l < m; ++l) {
					float coefficient1 = Q[iu] * C[w] * C[l];
					float integral     = 0.0;

#ifdef DEBUG
					std::cout << "Inner Coefficient: " << coefficient1 << std::endl;
#endif

					// Now we need to set up and calculate the integral
					// for the electron - nucleus interation.
					float A1 = A[w][0] + A[l][0];
					float A2 = A[w][1] + A[l][1];
					float A3 = A[w][2] + A[l][2];

					float s1 = s[w][0] + s[l][0];
					float s2 = s[w][1] + s[l][1];
					float s3 = s[w][2] + s[l][2];

					ElectronNucleiIntegrator integrator(this->size);
					integrator.A1 = A1;
					integrator.A2 = A2;
					integrator.A3 = A3;
					integrator.s1 = s1;
					integrator.s2 = s2;
					integrator.s3 = s3;
					integrator.R  = R[iu];

					integral = (float)integrator.Integrate(this->relerr, err);
					sum += integral * coefficient1;


#ifdef DEBUG
					std::cout << "Integral: " << integral << std::endl;
#endif
				}
			}
		}

		return -coefficient * sum;
	}

};