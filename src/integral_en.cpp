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
class ElectronNucleiIntegrator {

public:

	double lowerBounds[3];
	double upperBounds[3];
	GaussianWFN wavefn;

	int   maxCalls;
	float width;
	float epsilon;

	// Initializes an integral calculator.
	//     fn     = The wavefunction that the computation is being carried
	//              out for.
	//     width  = The number of half width half maxes to 
	//              limit the integral range to. Monte Carlo
	//              integration does not work for bounds at infinity.
	//
	//     ncalls = The maximum number of calls to the integrand function
	//              that the integrator is permitted to make in order to 
	//              calculate the integral. Higher values are slower and
	//              more accurate.
	ElectronNucleiIntegrator(GaussianWFN fn, float width, int ncalls) {
		lowerBounds    = new double[3];
		upperBounds    = new double[3];
		this->width    = width;
		this->maxCalls = ncalls;
		this->epsilon  = 0.25;
	}

	~ElectronNucleiIntegrator() {
		delete[] lowerBounds;
		delete[] upperBounds;
	}

	// Given the three indices that define the integral being calculated,
	// sets the lowerBounds and upperBounds members of the class to 
	// values appropriate for the given width, nuclei location, gaussian
	// center and gaussian width.
	void setBounds(int iu, int w, int l) {
		
	}

	// Integrates with the given relative error tolerance.
	float Integrate() {

		double result;
		double error;
		double lower[] = {-3/sqrtf(A1), -3/sqrtf(A1), -3/sqrtf(A1)};
		double upper[] = { 3/sqrtf(A1),  3/sqrtf(A1),  3/sqrtf(A1)};

		gsl_monte_function F = { 
			&ElectronNucleiIntegrator::integrand, 3, this
		};

		const gsl_rng_type *T;
  		gsl_rng *rng;

  		size_t calls = 3000;

		gsl_rng_env_setup();

		T   = gsl_rng_default;
		rng = gsl_rng_alloc(T);
 
		gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(3);

    	gsl_monte_vegas_integrate(&F, lower, upper, 3, calls, rng, s,
                               &result, &error);

		*err = error;

		return (float)result;
	}

private:
	

	static double integrand(double * x, size_t dim, void * params) {
		ElectronNucleiIntegrator * self = (ElectronNucleiIntegrator *)params;

		float x1 = x[0];
		float x2 = x[1];
		float x3 = x[2];

		float A1  = self->A1;
		float A2  = self->A2;
		float A3  = self->A3;
		float s1  = self->s1;
		float s2  = self->s2;
		float s3  = self->s3;
		float * R = self->R;

		float numerator = -(x1*x1)*A1 + x1*s1 - (x2*x2)*A2 + x2*s2 - (x3*x3)*A3 + x3*s3;
		numerator = expf(numerator);

		float d1 = R[0] - x1;
		float d2 = R[1] - x2;
		float d3 = R[2] - x3;

		float denominator = d1*d1 + d2*d2 + d3*d3;
		denominator = sqrtf(denominator);

		return numerator / denominator;
	}

};