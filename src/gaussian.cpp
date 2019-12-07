// Author:  Adam J. Robinson
// Summary: This file contains implementation of the Gaussian trial 
//          wavefunction used by MGSC. Some of the code is delegated to other
//          files for readability. Wherever possible, I attempt to name 
//          variables the same way that they are named in the mathematical
//          documentation.

#include <iostream>   // Used only in DEBUG mode for error messages
#include <mathimf.h>
#include <chrono>     // Used in profile mode 
#include "constants.h"
#include "integral_en.h"

// Stores all variables related to a Gaussian trial wavefunction and
// provides methods for calculating expectation values.
class GaussianWFN {
public:
	// These variables need to be assigned a value before the 
	// trial wavefunction can be used.
	int      m;  // Number of terms
	int      n;  // Number of inputs (3 * # of electrons)
	float ** A;  // Shape matrices (really vectors due to simplifications)
	float ** s;  // Shift matrices (really vectors)
	float  * C;  // Coefficients on each Gaussian term
	float ** R;  // List of nuclear coordinates
	int      Nu; // Number of nuclei
	float  * Q;  // List of nuclear charges

	// Initializes a Gaussian trial wavefunction with a given number of 
	// Gaussian terms. The variables annotated above need to be given a value
	// before the wavefunction is actually useable.
	GaussianWavefunction(int nterms, int nelectrons) {
		// In a future version, this will be 3 times and number of electrons.
		// Right now it is just going to be 3, corresponding to a single 
		// electron.
		n = nelectrons;
		m = nterms;
	}

public:

	// This section contains the primary functionality of the class,
	// calculation of the expectation value of the Hamiltonian.

	// Gets the expectation value of the Hamiltonian based on the
	// current set of parameters defining the wavefunction.
	float getHamiltonianExpectation() {
		float B = getNormalizationConstant();
		float T = getKineticExpectation(B);
		float V = getPotentialExpectation(B);

		return T + V;
	}

	// This calculates the value denoted 'B' in the documentation. 
	// The value is calculated such that |B|^2 multiplied by the 
	// integral over all space of |wavefunction|^2 will always be
	// unity.
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
				termsum  *= C[tidx] * C[outer_tidx];
				totalsum += termsum;
			}
		}		

		return 1.0 / sqrtf(totalsum);
	}

	// Given the appropriate normalization constant, returns the
	// expectation value of the kinetic energy for the given set 
	// of parameters.
	float getKineticExpectation(float B) {
		float coefficient  = -((B*B)*(N_hbar*N_hbar) / (2*N_me));
		      coefficient *= powrf(pi, ((float)n) / 2.0);

		float sum         = 0.0;

		for (int j_prime = 0; j_prime < m; ++j_prime) {
			for (int j = 0; j < m; ++j) {
				for (int i = 0; i < n; ++i) {
					float D_ijjp       = A[j_prime][i] + A[j][i];
					float F_ijjp       = s[j_prime][i] + s[j][i];
					float coefficient1 = C[j_prime] * C[j];

					float Gamma_t2       = (s[j][i]*s[j][i] - 2*A[j][i]);
					      Gamma_t2      *= (1/sqrtf(D_ijjp));
					float intermediate1  = (A[j][i]*A[j][i]);
					      intermediate1 *= (2*D_ijjp + F_ijjp*F_ijjp);
					      intermediate1 *= (1 / powrf(D_ijjp, 2.5));
					      Gamma_t2      += intermediate1;
					float intermediate2  = 2*A[j][i]*s[j][i]*F_ijjp;
						  intermediate2 *= (1 / powrf(D_ijjp, 1.5));
					      Gamma_t2      -= intermediate2;
					float Gamma          = expf((F_ijjp*F_ijjp) / (4*D_ijjp));
					      Gamma         *= Gamma_t2;

					// Now we compute the product term.
					float product = 1.0;
					for (int l = 0; l < n; ++l) {
						if (l != i) {
							float D_ljjp       = A[j_prime][l] + A[j][l];
							float F_ljjp       = s[j_prime][l] + s[j][l];
							product *= (1 / sqrtf(D_ljjp))
							product *= expf((F_ljjp*F_ljjp) / (4 * D_ljjp));
						}
					}

					sum += coefficient1 * Gamma * product;
				}
			}
		}

		return coefficient * sum;
	}

	

	// Given the appropriate normalization constant, returns the
	// expectation value of the potential energy for the given set 
	// of parameters.
	float getPotentialExpectation(float B) {
		float coefficient = (N_qe*B*B) / N_fpien;


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
				}
			}
		}

		return -coefficient * sum;
	}

}