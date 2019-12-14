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
#include "gaussian.h"



// Initializes a Gaussian trial wavefunction with a given number of 
// Gaussian terms. The variables annotated above need to be given a value
// before the wavefunction is actually useable.
GaussianWFN::GaussianWFN(int nterms, int nelectrons) {
	// In a future version, this will be 3 times and number of electrons.
	// Right now it is just going to be 3, corresponding to a single 
	// electron.
	n = nelectrons * 3;
	m = nterms;

	integrator = new ENIntegrator(this, width, ncalls);
}

GaussianWFN::~GaussianWFN() {
	delete integrator;
}
// This section contains the primary functionality of the class,
// calculation of the expectation value of the Hamiltonian.

// Gets the expectation value of the Hamiltonian based on the
// current set of parameters defining the wavefunction.
double GaussianWFN::getHamiltonianExpectation() {
	double B = getNormalizationConstant();
	double T = getKineticExpectation(B);
	double V = getPotentialExpectation(B);

	return T + V;
}

// This calculates the value denoted 'B' in the documentation. 
// The value is calculated such that |B|^2 multiplied by the 
// integral over all space of |wavefunction|^2 will always be
// unity.
double GaussianWFN::getNormalizationConstant() {
	// We basically need to sum the integral over all space of each
	// term. Each term is multiplied by it's corresponding C coefficient.
	// We then return the reciprocal of this sum.

	double totalsum = 0.0;
	for (int outer_tidx = 0; outer_tidx < m; ++outer_tidx) {
		for (int tidx = 0; tidx < m; ++tidx) {
			double termsum = 1.0;
			for (int iidx = 0; iidx < n; ++iidx) {
				double _s      = s[tidx][iidx] + s[outer_tidx][iidx];
				double _A      = A[tidx][iidx] + A[outer_tidx][iidx];
				double ssquare = _s*_s;
				double aterm   = _A;
				termsum *= sqrt(pi / aterm)*exp(ssquare / (4.0 * aterm));
			}
			termsum  *= C[tidx] * C[outer_tidx];
			totalsum += termsum;
		}
	}		

	return 1.0 / sqrt(totalsum);
}

// Given the appropriate normalization constant, returns the
// expectation value of the kinetic energy for the given set 
// of parameters.
double GaussianWFN::getKineticExpectation(double B) {
	double coefficient  = -((B*B)*(N_hbar*N_hbar) / (2*N_me));
	       coefficient *= pow(pi, ((double)n) / 2.0);

	double sum         = 0.0;

	for (int j_prime = 0; j_prime < m; ++j_prime) {
		for (int j = 0; j < m; ++j) {
			for (int i = 0; i < n; ++i) {
				double D_ijjp       = A[j_prime][i] + A[j][i];
				double F_ijjp       = s[j_prime][i] + s[j][i];
				double coefficient1 = C[j_prime] * C[j];

				double Gamma_t2       = (s[j][i]*s[j][i] - 2*A[j][i]);
				       Gamma_t2      *= (1/sqrt(D_ijjp));
				double intermediate1  = (A[j][i]*A[j][i]);
				       intermediate1 *= (2*D_ijjp + F_ijjp*F_ijjp);
				       intermediate1 *= (1 / pow(D_ijjp, 2.5));
				       Gamma_t2      += intermediate1;
				double intermediate2  = 2*A[j][i]*s[j][i]*F_ijjp;
					   intermediate2 *= (1 / pow(D_ijjp, 1.5));
				       Gamma_t2      -= intermediate2;
				double Gamma          = exp((F_ijjp*F_ijjp) / (4*D_ijjp));
				       Gamma         *= Gamma_t2;

				// Now we compute the product term.
				double product = 1.0;
				for (int l = 0; l < n; ++l) {
					if (l != i) {
						double D_ljjp       = A[j_prime][l] + A[j][l];
						double F_ljjp       = s[j_prime][l] + s[j][l];
						product *= (1 / sqrt(D_ljjp));
						product *= exp((F_ljjp*F_ljjp) / (4 * D_ljjp));
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
double GaussianWFN::getPotentialExpectation(double B) {
	double coefficient = (N_qe*B*B) / N_fpien;

	double sum = 0.0;
	for (int iu = 0; iu < Nu; ++iu) {
		// Normally there would be another loop here over 
		// the number of electrons, but I'm omitting that in this
		// version because of time constraints. I am only doing
		// one electron systems for now.

		for (int w = 0; w < m; ++w) {
			for (int l = 0; l < m; ++l) {
				double coefficient1 = Q[iu] * C[w] * C[l];
				double integral     = 0.0;

				integral = integrator->Integrate(iu, w, l);
				sum += integral * coefficient1;
			}
		}
	}

	return -coefficient * sum;
}