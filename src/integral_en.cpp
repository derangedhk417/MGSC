// Author:  Adam J. Robinson
// Summary: This file contains implementation of the integral that 
//          characterizes the interaction between nuclei and electrons,
//          as defined in the mathematical documentation. It used the
//          VEGAS Monte Carlo integration algorithm provided by GSL to
//          accomplish this task efficiently.

#include <iostream>   // Used only in DEBUG mode for error messages
#include <mathimf.h>
#include <chrono>     // Used in profile mode 
#include "constants.h"
#include "gaussian.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_errno.h>

// Performs the Electron Nuclei integral for the given 
// Gaussian wavefunction.
class ENIntegrator {

public:

	gsl_monte_vegas_state * state;

	double lowerBounds[3];
	double upperBounds[3];
	GaussianWFN * wavefn;

	int    maxCalls;
	double width;
	double epsilon;

	// Initializes an integral calculator.
	//     fn     = The wavefunction that the computation is being carried
	//              out for.
	//     width  = The number of half width half maxes to 
	//              limit the integral range to. Monte Carlo
	//              integration does not work for bounds at infinity.
	//              typically, 3 < width < 6
	//
	//     ncalls = The maximum number of calls to the integrand function
	//              that the integrator is permitted to make in order to 
	//              calculate the integral. Higher values are slower and
	//              more accurate.
	ENIntegrator(GaussianWFN * fn, double width, int ncalls) {
		lowerBounds    = new double[3];
		upperBounds    = new double[3];
		this->wavefn   = fn;
		this->width    = width;
		this->maxCalls = ncalls;
		this->epsilon  = 0.25;

		state = gsl_monte_vegas_alloc(3);
	}

	~ENIntegrator() {
		delete[] lowerBounds;
		delete[] upperBounds;
		gsl_monte_vegas_free(state);
	}

	// These are set temporarily before an integral is performed.
	double A1;
	double A2;
	double A3;

	double s1;
	double s2;
	double s3;

	double R1;
	double R2;
	double R3;

	// Given the three indices that define the integral being calculated,
	// sets the lowerBounds and upperBounds members of the class to 
	// values appropriate for the given width, nuclei location, gaussian
	// center and gaussian width.
	void setBounds(int iu, int w, int l) {
		A1 = wavefn->A[w][0] + wavefn->A[l][0];
		A2 = wavefn->A[w][1] + wavefn->A[l][1];
		A3 = wavefn->A[w][2] + wavefn->A[l][2];

		s1 = wavefn->s[w][0] + wavefn->s[l][0];
		s2 = wavefn->s[w][1] + wavefn->s[l][1];
		s3 = wavefn->s[w][2] + wavefn->s[l][2];

		R1 = wavefn->R[iu][0];
		R2 = wavefn->R[iu][1];
		R3 = wavefn->R[iu][2];

		// Here we calculate the appropriate bounds.
		double bound_xrl = R1 - epsilon;
		double bound_xel = - (s1 / (2*A1)) - width*sqrt(ln2 / A1);

		double bound_xru = R1 + epsilon;
		double bound_xeu = - (s1 / (2*A1)) + width*sqrt(ln2 / A1);

		lowerBounds[0] = fmin(bound_xrl, bound_xel);
		upperBounds[0] = fmax(bound_xru, bound_xeu);

		bound_xrl = R2 - epsilon;
		bound_xel = - (s2 / (2*A2)) - width*sqrt(ln2 / A2);

		bound_xru = R2 + epsilon;
		bound_xeu = - (s2 / (2*A2)) + width*sqrt(ln2 / A2);

		lowerBounds[1] = fmin(bound_xrl, bound_xel);
		upperBounds[1] = fmax(bound_xru, bound_xeu);

		bound_xrl = R3 - epsilon;
		bound_xel = - (s3 / (2*A3)) - width*sqrt(ln2 / A3);

		bound_xru = R3 + epsilon;
		bound_xeu = - (s3 / (2*A3)) + width*sqrt(ln2 / A3);

		lowerBounds[2] = fmin(bound_xrl, bound_xel);
		upperBounds[2] = fmax(bound_xru, bound_xeu);
	}

	double Integrate(int iu, int w, int l) {
		setBounds(iu, w, l);

		double result;
		double error;

		gsl_monte_function F = { 
			&ENIntegrator::integrand, 3, this
		};

		const gsl_rng_type *T;
  		gsl_rng *rng;

  		size_t calls = (size_t)maxCalls;

		gsl_rng_env_setup();

		T   = gsl_rng_default;
		rng = gsl_rng_alloc(T);
 
 		// Reinitialize the state so that the data from the previous
 		// integration will be removed.
 		gsl_monte_vegas_init(state);

    	gsl_monte_vegas_integrate(
    		&F, lowerBounds, upperBounds, 3, calls, 
    		rng, state, &result, &error
		);

		return result;
	}

private:

	static double integrand(double * x, size_t dim, void * params) {
		ENIntegrator * self = (ENIntegrator *)params;

		float x1 = x[0];
		float x2 = x[1];
		float x3 = x[2];

		float A1  = self->A1;
		float A2  = self->A2;
		float A3  = self->A3;
		float s1  = self->s1;
		float s2  = self->s2;
		float s3  = self->s3;
		float R1  = self->R1;
		float R2  = self->R2;
		float R3  = self->R3;

		float numerator  = -(x1*x1)*A1 + x1*s1 - (x2*x2)*A2 + x2*s2 
		      numerator += -(x3*x3)*A3 + x3*s3;
		numerator = expf(numerator);

		float d1 = R1 - x1;
		float d2 = R2 - x2;
		float d3 = R3 - x3;

		float denominator = d1*d1 + d2*d2 + d3*d3;
		      denominator = sqrtf(denominator);

		return numerator / denominator;
	}

};